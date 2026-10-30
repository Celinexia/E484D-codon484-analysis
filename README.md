# Codon 484 analysis for SARS-CoV-2 E484D: Figure 10, panels A and B

Scripts and tables for Figure 10A (codon 484 substitution landscape) and 10B (E484D allele fraction per SRA run) of the TMEM106B / SARS-CoV-2 SE484D manuscript. Raw sequencing reads are not included. Accession numbers are in the manuscript Data Availability section.

## Panel A

- **Script:** `figure10_panelA_color.R`. Run `Rscript figure10_panelA_color.R`. Writes `Figure10-panelA-color.pdf` and `.png` (600 dpi).
- **Input:** `lapis_codon484_counts.csv`. The counts are typed into the script, not read from this file. The file records the same values for reference.
- **Source:** LAPIS/CoV-Spectrum, GenBank INSDC genomes, 12 Oct 2026. The query itself is not included. Values were transcribed from the manuscript Methods.
- **Colours:** Omicron, Beta/Gamma/Mu and Kappa are coloured by lineage. E484D is teal (`#0F9D9A`), as in the legend. Blue and orange are reserved for the 23014 A>C and A>T routes in panel B.
- **Software:** R 4.3.1, ggplot2 3.4.4, scales 1.2.1.

## Panel B

- **Script:** `figure10_panelB_regrouped.R`. Run `Rscript figure10_panelB_regrouped.R`. Writes `Figure10-panelB-regrouped.pdf` and `.png`.
- **Input:** `panelB_data_regrouped.csv`, 406 rows, one per SRA run. The MA1 mouse-adaptation points (0.96 lung, 0.92 supernatant) are typed into the script, from Yan et al.
- **Layout:** five columns, IN VIVO / EX VIVO / NOT RECORDED. Helvetica font, metric-compatible with Arial.
- **Software:** R 4.3.1, ggplot2 3.4.4, scales 1.2.1.

### How the panel B inputs are built

1. `fetch_sra_metadata.py` queries NCBI E-utilities for the metadata of the runs in `ncbi_sra_E484D_matches.csv`, and writes `sra_run_metadata.csv`. Set `NCBI_EMAIL` first.
2. `plot_figure10_final.py` classifies each run into a context (its `context()` function) and prints the counts. Its output matches `panelB_data.csv` exactly: 0 of 406 runs differ. The script also draws a full draft figure, which is not the submitted layout.
3. `make_panelB_regrouped.py` rebuilds `panelB_data_regrouped.csv` from `panelB_data.csv` and `sra_run_metadata.csv`. It changes only the Treatment-selection runs:
   - PRJEB55794 (the monoclonal-antibody run) goes to Clinical surveillance, as the Methods state.
   - PRJNA1035781 and PRJNA692078 (drug-selection passages) go to Drug selection (cell culture).
   The rebuilt file matches the shipped one, except for floating-point noise in the 16th digit of `vaf`.

### Context counts (manuscript)

Clinical surveillance 227 (226 + the mAb run). Animal 17 (15 experimental + 1 wild deer + the MA1 lung point). Cell culture and drug selection 39 (15 + 23 + the MA1 supernatant point). Wastewater 14. Unspecified 111. Total 408 = 406 runs + 2 MA1 points.

## Run order

```bash
export NCBI_EMAIL=you@example.org
python3 fetch_sra_metadata.py                 # optional: rebuilds sra_run_metadata.csv
python3 plot_figure10_final.py                # optional: prints context counts
python3 make_panelB_regrouped.py              # rebuilds panelB_data_regrouped.csv
Rscript figure10_panelA_color.R
Rscript figure10_panelB_regrouped.R
```

## Known gaps

- The LAPIS/CoV-Spectrum query behind panel A is not in this folder.
- Run SRR18915556 has no metadata in NCBI SRA or ENA, so it is classed as Unspecified by default. The manuscript counts it.
