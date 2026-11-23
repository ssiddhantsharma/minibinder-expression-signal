"""How much signal is in display-sort enrichment? Almost none -- for anything.

Run from the repo root:
  uvx --python 3.12 --with matplotlib,pandas,numpy,scipy python scripts/plot_fimh_label_signal.py
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

B, rng = 2000, np.random.default_rng(0)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "STIX Two Text", "DejaVu Serif"],
    "axes.labelsize": 14, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "axes.linewidth": 1.0, "axes.spines.top": False, "axes.spines.right": False,
})
BLUE, SLATE, DARK = "#1f6aa5", "#8d99a6", "#111417"
KD = dict(zip("AVLIPFMWGSTCYNQDEKRH",
              [1.8, 4.2, 3.8, 4.5, -1.6, 2.8, 1.9, -.9, -.4, -.8,
               -.7, 2.5, -1.3, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5, -3.2]))

p = pd.read_csv("data/fimh_monoyield_predictions.csv")
t = pd.read_csv("data/fimh_display_expression.csv")
t["name"] = "fimh_" + t.index.astype(str)
d = p.merge(t[["name", "sequence", "log2_expressed", "expression_flag"]], on="name",
            validate="one_to_one")
# the join key is positional, so check it against what was actually scored
assert (d.raw_aa_seq == d.sequence).all(), "prediction/expression rows are misaligned"
ok = d[d.expression_flag == "ok"].copy()
s = ok.raw_aa_seq

feats = {
    "net charge  (K+R − D+E)": s.map(lambda x: sum(c in "KR" for c in x) - sum(c in "DE" for c in x)),
    "aromatic fraction": s.map(lambda x: sum(c in "FWY" for c in x) / len(x)),
    "hydrophobicity  (GRAVY)": s.map(lambda x: np.mean([KD.get(c, 0) for c in x])),
    "ESMC-MonoYield  (6B parameters)": ok.monoyield_log10_mgL,
    "ILV fraction": s.map(lambda x: sum(c in "ILV" for c in x) / len(x)),
}
y = ok.log2_expressed.values
idx = rng.integers(0, len(y), (B, len(y)))


def rho2_ci(v):
    v = np.asarray(v)
    bs = np.array([spearmanr(v[i], y[i]).statistic for i in idx])
    return 100 * spearmanr(v, y).statistic ** 2, *np.percentile(100 * bs ** 2, [2.5, 97.5])


rows = sorted(((k, *rho2_ci(v)) for k, v in feats.items()), key=lambda r: -r[1])

fig, ax = plt.subplots(figsize=(10.5, 6.2))
names = [r[0] for r in rows]
r2 = [r[1] for r in rows]
err = np.array([[r[1] - r[2] for r in rows], [r[3] - r[1] for r in rows]])
cols = [BLUE if "MonoYield" in n else SLATE for n in names]
pos = np.arange(len(names))
ax.barh(pos, r2, 0.64, color=cols)
ax.errorbar(r2, pos, xerr=err, fmt="none", ecolor=DARK, elinewidth=1.3, capsize=4)
for i, r in enumerate(rows):
    ax.text(r[3] + 0.06, i, f"{r[1]:.1f}%", va="center", fontsize=14,
            color=BLUE if "MonoYield" in names[i] else DARK,
            fontweight="bold" if "MonoYield" in names[i] else "normal")
ax.set_yticks(pos)
ax.set_yticklabels(names, fontsize=13.5)
for lab, n in zip(ax.get_yticklabels(), names):
    if "MonoYield" in n:
        lab.set_color(BLUE); lab.set_fontweight("bold")
ax.invert_yaxis()
ax.set_xlim(0, 2.0)
ax.set_xlabel("rank variance in measured expression explained  "
              "(Spearman $\\rho^2$, 95% bootstrap CI)", labelpad=10)

ax.set_title("Nothing predicts display-sort enrichment \u2014 not a 6-billion-parameter\n"
             "yield model, not any simple sequence feature",
             fontsize=16, pad=16, color=DARK)
fig.text(0.5, 0.075,
         "5,114 de novo minibinders, all 85 aa — Chazin-Gray et al., Nat Commun 2026 · "
         "scored with ESMC-MonoYield (Setsu et al.)",
         ha="center", fontsize=11.5, color=DARK)
fig.text(0.5, 0.005,
         "MonoYield's interval overlaps every hand-coded feature's \u2014 it is not distinguishable from any of them.\n"
         "Enrichment from cDNA/yeast display, not the E. coli mg/L MonoYield predicts: the label is what fails, not the model.",
         ha="center", va="bottom", fontsize=11, color=DARK, linespacing=1.5)
fig.tight_layout(rect=[0, 0.115, 1, 1])
fig.savefig("figures/fimh_label_signal.png", dpi=220)
print("wrote figures/fimh_label_signal.png")
