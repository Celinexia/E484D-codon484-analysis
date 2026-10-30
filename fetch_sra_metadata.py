"""Fetch BioProject / study metadata for each E484D run from NCBI E-utilities.

Writes sra_run_metadata.csv next to the input. Cached: re-running skips the
fetch if the output already covers every run.
"""

import os
import subprocess
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET

import pandas as pd

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 80
IN_CSV = "ncbi_sra_E484D_matches.csv"
OUT_CSV = "sra_run_metadata.csv"


def fetch(accessions):
    params = {"db": "sra", "id": ",".join(accessions), "rettype": "xml",
              "tool": "e484d-manuscript", "email": os.environ.get("NCBI_EMAIL", "")}
    url = f"{EUTILS}?{urllib.parse.urlencode(params)}"
    # curl, not urllib: this network terminates TLS with its own CA, which the
    # system keychain trusts but Python's bundled certifi does not
    for attempt in range(4):
        try:
            out = subprocess.run(
                ["curl", "-sS", "--fail", "--max-time", "120", url],
                capture_output=True, check=True).stdout
            return ET.fromstring(out)
        except Exception as exc:  # transient 429/500 from NCBI are common
            if attempt == 3:
                raise
            print(f"  retry {attempt + 1} after {exc}", file=sys.stderr)
            time.sleep(3 * (attempt + 1))


def text(node, path):
    el = node.find(path) if node is not None else None
    return el.text.strip() if el is not None and el.text else None


def parse(root):
    rows = []
    for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
        study = pkg.find("STUDY")
        bioproject = None
        if study is not None:
            for eid in study.findall(".//EXTERNAL_ID"):
                if eid.get("namespace") == "BioProject":
                    bioproject = eid.text
                    break
        common = {
            "bioproject": bioproject,
            "srp": (study.get("accession") if study is not None else None),
            "study_title": text(study, ".//STUDY_TITLE"),
            "study_abstract": text(study, ".//STUDY_ABSTRACT"),
            "library_strategy": text(pkg, ".//LIBRARY_STRATEGY"),
            "library_source": text(pkg, ".//LIBRARY_SOURCE"),
            "biosample": next(
                (e.text for e in pkg.findall(".//EXTERNAL_ID")
                 if e.get("namespace") == "BioSample"), None),
            "sample_title": text(pkg, ".//SAMPLE/TITLE"),
            "scientific_name": text(pkg, ".//SAMPLE_NAME/SCIENTIFIC_NAME"),
            "center": (pkg.find("SUBMISSION").get("center_name")
                       if pkg.find("SUBMISSION") is not None else None),
        }
        attrs = {}
        for a in pkg.findall(".//SAMPLE_ATTRIBUTE"):
            tag, val = text(a, "TAG"), text(a, "VALUE")
            if tag:
                attrs[tag.lower()] = val
        for key in ("isolation_source", "isolate", "passage_history",
                    "sample_type", "tissue", "cell_line", "source_name"):
            common[key] = attrs.get(key)
        for run in pkg.findall(".//RUN"):
            acc = run.get("accession")
            if acc:
                rows.append({"Run": acc, **common})
    return rows


def main():
    runs = pd.read_csv(IN_CSV).Run.dropna().unique().tolist()
    try:
        have = pd.read_csv(OUT_CSV)
        if set(runs) <= set(have.Run):
            print(f"cached: {OUT_CSV} already covers {len(runs)} runs")
            return
    except FileNotFoundError:
        pass

    rows = []
    for i in range(0, len(runs), BATCH):
        chunk = runs[i:i + BATCH]
        print(f"fetching {i + 1}-{i + len(chunk)} of {len(runs)}")
        rows += parse(fetch(chunk))
        time.sleep(0.5)  # stay under the 3 req/s anonymous limit

    out = pd.DataFrame(rows).drop_duplicates(subset="Run")
    out = out[out.Run.isin(runs)]
    out.to_csv(OUT_CSV, index=False)
    missing = sorted(set(runs) - set(out.Run))
    print(f"wrote {OUT_CSV}: {len(out)}/{len(runs)} runs, "
          f"{out.bioproject.nunique()} BioProjects")
    if missing:
        print(f"no metadata returned for {len(missing)}: {missing[:10]}")


if __name__ == "__main__":
    main()
