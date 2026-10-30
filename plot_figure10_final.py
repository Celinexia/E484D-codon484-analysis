"""Figure 10. In vivo evolution at spike residue 484.

Answers Reviewer #1 point 2 (E484D does not replicate in ACE2-deficient mice,
so show in vivo relevance) by separating two claims the critique conflates:
whether TMEM106B supports *spread* in an ACE2-null host, and whether E484D is
an accessible in vivo evolutionary *outcome*. This figure addresses the second.

  A  codon 484 (GAA) reaches six substitutions in one nucleotide step. E484D is
     the only one reachable by two independent routes, yet is the rarest in
     nature -- so its rarity is selective, not mutational.
  B  every setting in which E484D has been observed, in vivo above the rule and
     ex vivo below it, on one common axis: the fraction of reads carrying the
     substitution. SRA runs and the MA1 mouse-adaptation lineage plot together
     because both are read fractions.

Colour discipline, inherited from plot_e484d.py: blue/orange mean the
nucleotide route at 23014 and nothing else. The MA1 points are a different
study, so they are keyed by marker shape, not by hue.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle

CSV = "ncbi_sra_E484D_matches.csv"
META = "sra_run_metadata.csv"
OUT = "Figure10-final"

BLUE = "#2a78d6"     # 23014 A>C
ORANGE = "#eb6834"   # 23014 A>T
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
FAINT = "#d8d7d0"
RULE = "#b9b8b0"
PANEL = "#f4f3ef"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.edgecolor": MUTED,
    "axes.linewidth": 0.6,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "text.color": INK,
    "axes.labelcolor": INK,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
})

N_GENOMES = 9_438_404
NEIGHBOURS = [
    # label, codon(s), base index, nt change, genomes in nature
    ("E484K", "AAA", 0, "23012 G>A", 370_685),
    ("E484Q", "CAA", 0, "23012 G>C", 8_237),
    ("E484A", "GCA", 1, "23013 A>C", 4_306_700),
    ("E484G", "GGA", 1, "23013 A>G", 170),
    ("E484V", "GTA", 1, "23013 A>T", None),
    ("E484D", "GAC / GAT", 2, "23014 A>C or A>T", 229),
]

# MA1 mouse-adaptation lineage, Yan et al. Virus Evolution 2022, Suppl.
# Dataset 1: percent of reads carrying E484D. Same quantity as the SRA allele
# fractions, so the two plot on one axis.
MA1 = {"Mouse lung, MA1 adaptation": 0.96,
       "Cell supernatant, MA1": 0.92}

# ---------------------------------------------------------------------------
TREATMENT = {"PRJNA1035781", "PRJNA692078", "PRJEB55794"}
LAB_PROJECTS = {"PRJNA804321", "PRJNA1018257"}


def context(r):
    src = str(r.get("isolation_source") or "").lower()
    title = str(r.get("study_title") or "").lower()
    host = str(r.get("host") or "")
    if r.get("SRA_organism") == "wastewater metagenome" or "wastewater" in title:
        return "Wastewater"
    if any(k in src for k in ("vero", "cell culture", "hek", "cell line",
                              "supernatant")) or host == "Chlorocebus aethiops":
        return "Cell culture"
    if host in ("Hamster", "Mus musculus") or any(
            k in title for k in ("hamster", "mouse adapted", "serial passag")):
        return "Experimental animal"
    if host == "Odocoileus virginianus" or "free-ranging" in title:
        return "Wild animal"
    if r.get("bioproject") in LAB_PROJECTS:
        return ("Experimental animal" if "mouse" in src or host == "Mus musculus"
                else "Cell culture")
    if r.get("bioproject") in TREATMENT:
        return "Treatment-selection"
    if host == "Homo sapiens" or any(k in src for k in (
            "swab", "clinical", "patient", "nasal", "nasopharyn", "sputum")):
        return "Clinical surveillance"
    return "Unspecified"


df = pd.read_csv(CSV).merge(pd.read_csv(META), on="Run", how="left")
df["vaf"] = df.alternate_read_depth / df.read_depth
df["context"] = df.apply(context, axis=1)

# row order, top to bottom, with the in vivo / ex vivo rule between them
ROWS = [
    ("group", "IN VIVO"),
    ("sra", "Clinical surveillance", "Clinical surveillance", "human"),
    ("sra", "Wild white-tailed deer", "Wild animal", "animal"),
    ("sra", "Experimental animal", "Experimental animal", "animal"),
    ("ma1", "Mouse lung, MA1 adaptation", None, "animal"),
    ("sra", "Treatment-selection", "Treatment-selection", "human"),
    ("group", "EX VIVO"),
    ("sra", "Cell culture", "Cell culture", None),
    ("ma1", "Cell supernatant, MA1", None, None),
    ("sra", "Wastewater", "Wastewater", None),
    ("group", "CONTEXT NOT RECORDED"),
    ("sra", "Unspecified", "Unspecified", None),
]

# ---------------------------------------------------------------------------
# Layout -- 180 mm two-column width
# ---------------------------------------------------------------------------
W = 7.09
A_H = 2.70
B_H = 2.72
BOT = 1.15
GAP = 0.70
TOP = 0.26
H = BOT + B_H + GAP + A_H + TOP

fig = plt.figure(figsize=(W, H))


def ax_in(x, y, w, h):
    return fig.add_axes([x / W, y / H, w / W, h / H])


ax_a = ax_in(0.30, BOT + B_H + GAP, W - 0.60, A_H)
ax_b = ax_in(1.62, BOT, 4.62, B_H)


def panel_tag(ax, letter, dx, dy):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(dx, dy), textcoords="offset points",
                fontsize=9, fontweight="bold", va="bottom", ha="left",
                color=INK, annotation_clip=False)


# ===========================================================================
# A -- mutational access vs selective retention
# ===========================================================================
ax_a.set_xlim(0, 100)
ax_a.set_ylim(0, 100)
ax_a.axis("off")
panel_tag(ax_a, "A", -16, 0)

ax_a.text(0, 104, "Codon 484 reaches six substitutions in one step; "
                  "which one is retained depends on the receptor context",
          fontsize=7.8, fontweight="bold", color=INK, va="top", ha="left")

cx, cy, bw, bh = 3.0, 56.0, 5.4, 10.5
BASE_POS = [23012, 23013, 23014]
for i, b in enumerate("GAA"):
    x = cx + i * bw
    hot = i == 2
    ax_a.add_patch(Rectangle((x, cy), bw, bh, facecolor="white",
                             edgecolor=INK if hot else MUTED,
                             linewidth=1.1 if hot else 0.7, zorder=3))
    ax_a.text(x + bw / 2, cy + bh / 2, b, fontsize=9,
              fontweight="bold" if hot else "normal",
              color=INK if hot else INK2, va="center", ha="center", zorder=4)
    ax_a.text(x + bw / 2, cy - 2.0, str(BASE_POS[i]), fontsize=5.6,
              color=MUTED, va="top", ha="center")
ax_a.text(cx + 1.5 * bw, cy + bh + 2.6, "E484  (Glu)", fontsize=7.2,
          fontweight="bold", color=INK, va="bottom", ha="center")
ax_a.text(cx + 1.5 * bw, cy - 9.5, "parental codon", fontsize=6.3,
          color=INK2, va="top", ha="center")

x_sub, x_cod, x_chg, x_cnt = 34.0, 47.0, 58.0, 79.0
rows_a = np.linspace(88, 36, len(NEIGHBOURS))
base_x = [cx + (i + 0.5) * bw for i in range(3)]

for xx, head in ((x_sub, "Substitution"), (x_cod, "Codon"),
                 (x_chg, "Nucleotide change"),
                 (x_cnt, f"Genomes in nature (of {N_GENOMES:,})")):
    ax_a.text(xx, 97, head, fontsize=5.9, color=INK2, va="center", ha="left")
ax_a.plot([x_sub - 1.5, 100], [93.5, 93.5], color=RULE, linewidth=0.6, zorder=1)

for (lab, codon, bi, change, cnt), yy in zip(NEIGHBOURS, rows_a):
    hi = lab == "E484D"
    col = INK if hi else INK2
    rads = (0.10, -0.24) if hi else (0.10,)
    for rad in rads:
        ax_a.add_patch(FancyArrowPatch(
            (base_x[bi], cy + bh if yy > cy else cy), (x_sub - 1.5, yy),
            connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
            mutation_scale=6, linewidth=1.15 if hi else 0.6,
            color=col if hi else FAINT, zorder=2, shrinkA=2, shrinkB=2))
    ax_a.text(x_sub, yy, lab, fontsize=7.3,
              fontweight="bold" if hi else "normal", color=col,
              va="center", ha="left")
    ax_a.text(x_cod, yy, codon, fontsize=6.4, color=MUTED,
              va="center", ha="left")
    ax_a.text(x_chg, yy, change, fontsize=6.0,
              fontweight="bold" if hi else "normal",
              color=col if hi else MUTED, va="center", ha="left")
    if cnt is None:
        txt = "not observed"
    else:
        pct = 100 * cnt / N_GENOMES
        p = (f"{pct:.4f}%" if pct < 0.01 else
             f"{pct:.3f}%" if pct < 1 else f"{pct:.1f}%")
        txt = f"{cnt:,}   {p}"
    ax_a.text(x_cnt, yy, txt, fontsize=6.2,
              fontweight="bold" if hi else "normal",
              color=col if hi else MUTED, va="center", ha="left")

ax_a.text(x_sub - 3.0, rows_a[-1] - 8.5, "two independent routes",
          fontsize=6.0, fontweight="bold", color=INK, va="center", ha="left")

ax_a.add_patch(Rectangle((0, 0), 100, 22.0, facecolor=PANEL,
                         edgecolor="none", zorder=1))
ax_a.text(2.0, 19.0, "E484D is the only substitution at this residue reachable "
                    "by two independent nucleotide changes, yet it is the "
                    "rarest observed.",
          fontsize=6.9, fontweight="bold", color=INK, va="top", ha="left")
ax_a.text(2.0, 13.5, "Its rarity is therefore selective, not mutational: where "
                    "ACE2 is engaged efficiently E484A/K/Q are retained and "
                    "define VOC lineages, while E484D\nis retained where that "
                    "constraint is relaxed — the condition under which "
                    "TMEM106B-dependent entry becomes available.",
          fontsize=6.5, color=INK2, va="top", ha="left", linespacing=1.45)

# ===========================================================================
# B -- every setting E484D has been observed in, in vivo above ex vivo
# ===========================================================================
rng = np.random.default_rng(11)
y_pos, y_lab, group_spans = {}, {}, []
y = 0.0
for item in ROWS:
    if item[0] == "group":
        y -= 0.72
        group_spans.append((item[1], y))
        y -= 0.30
    else:
        y -= 1.0
        y_pos[item[1]] = y
        y_lab[item[1]] = item

y_bottom = y - 0.75
ax_b.set_ylim(y_bottom, 0.25)
ax_b.set_xlim(-0.02, 1.16)
ax_b.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax_b.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.0"])
ax_b.set_xlabel("Fraction of reads carrying E484D")
ax_b.set_yticks([])
ax_b.spines[["top", "right", "left"]].set_visible(False)
ax_b.tick_params(length=3, width=0.6, color=MUTED)
ax_b.set_xbound(-0.02, 1.16)

# consensus reference
ax_b.axvline(0.8, color=RULE, linewidth=0.7, linestyle=(0, (4, 3)), zorder=1)
ax_b.text(0.8, 0.30, "consensus level", fontsize=6.0, color=INK2,
          va="bottom", ha="center")

# group headers and the in vivo / ex vivo rule
for name, yy in group_spans:
    ax_b.text(-0.055, yy, name, fontsize=6.3, fontweight="bold", color=INK,
              va="center", ha="left", clip_on=False)
    ax_b.plot([-0.02, 1.08], [yy - 0.34, yy - 0.34], color=RULE,
              linewidth=0.6, zorder=1, clip_on=False)

for label, (kind, lab, ctx, _h) in y_lab.items():
    yy = y_pos[label]
    ax_b.axhline(yy, color="#efeee9", linewidth=0.5, zorder=0)
    ax_b.text(-0.055, yy, lab, fontsize=6.5, color=INK2, va="center",
              ha="right", clip_on=False)

    if kind == "sra":
        g = df[df.context == ctx]
        for allele, col in (("C", BLUE), ("T", ORANGE)):
            v = g.loc[g.alternate == allele, "vaf"].values
            if not len(v):
                continue
            jit = rng.uniform(-0.19, 0.19, len(v))
            ax_b.scatter(v, yy + jit, s=7, facecolor=col, edgecolor="none",
                         alpha=0.55, zorder=3, rasterized=True)
        med = g.vaf.median()
        ax_b.plot([med, med], [yy - 0.30, yy + 0.30], color=INK,
                  linewidth=1.3, zorder=5, solid_capstyle="butt")
        ax_b.text(1.10, yy, f"{len(g)}", fontsize=6.2, color=INK2,
                  va="center", ha="right")
    else:
        v = MA1[lab]
        ax_b.scatter([v], [yy], s=30, marker="D", facecolor="white",
                     edgecolor=INK, linewidth=1.0, zorder=5)
        ax_b.text(1.10, yy, "1", fontsize=6.2, color=INK2,
                  va="center", ha="right")

ax_b.text(1.10, 0.30, "n", fontsize=6.2, color=INK2, va="bottom", ha="right")

ax_b.set_title("E484D arises in vivo and ex vivo, and reaches consensus in "
               "natural human infection",
               fontsize=7.8, fontweight="bold", color=INK, pad=16,
               loc="left", x=-0.265)
panel_tag(ax_b, "B", -104, 20)

handles = [
    Line2D([], [], marker="o", linestyle="none", markersize=3.4,
           markerfacecolor=BLUE, markeredgecolor="none",
           label="23014 A>C   (gaA/gaC)"),
    Line2D([], [], marker="o", linestyle="none", markersize=3.4,
           markerfacecolor=ORANGE, markeredgecolor="none",
           label="23014 A>T   (gaA/gaT)"),
    Line2D([], [], color=INK, linewidth=1.3, label="median"),
    Line2D([], [], marker="D", linestyle="none", markersize=4,
           markerfacecolor="white", markeredgecolor=INK, markeredgewidth=1.0,
           label="MA1 lineage (Yan et al. 2022)"),
]
ax_b.legend(handles=handles, loc="upper left",
            bbox_to_anchor=(-0.265, -0.165), frameon=False, ncol=4,
            handlelength=1.1, handletextpad=0.45, columnspacing=1.5,
            fontsize=6.2, borderpad=0)

ax_b.annotate(
    "SRA runs carrying E484D at genome position 23014 (n = 406), grouped by "
    "how the sample arose; one point per run, horizontally jittered.\n"
    "MA1 points are a separate study — percent of reads carrying E484D in a "
    "mouse-adapted lineage and its cell-culture passage.",
    xy=(-0.265, -0.300), xycoords="axes fraction", fontsize=5.9, color=MUTED,
    va="top", ha="left", linespacing=1.45, annotation_clip=False)

fig.savefig(f"{OUT}.pdf", dpi=600)
fig.savefig(f"{OUT}.png", dpi=300)
print(f"wrote {OUT}.pdf / {OUT}.png  ({W:.2f} x {H:.2f} in)")
for item in ROWS:
    if item[0] == "sra":
        g = df[df.context == item[2]]
        print(f"  {item[1]:28s} n={len(g):3d}  median={g.vaf.median():.3f}")
