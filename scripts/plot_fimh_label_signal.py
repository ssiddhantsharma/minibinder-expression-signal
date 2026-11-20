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
d = p.merge(t[["name", "log2_expressed", "expression_flag"]], on="name", validate="one_to_one")
ok = d[d.expression_flag == "ok"].copy()
s = ok.raw_aa_seq

feats = {
    "net charge  (K+R − D+E)": s.map(lambda x: sum(c in "KR" for c in x) - sum(c in "DE" for c in x)),
    "aromatic fraction": s.map(lambda x: sum(c in "FWY" for c in x) / len(x)),
    "hydrophobicity  (GRAVY)": s.map(lambda x: np.mean([KD.get(c, 0) for c in x])),
    "ESMC-MonoYield  (6B parameters)": ok.monoyield_log10_mgL,
    "ILV fraction": s.map(lambda x: sum(c in "ILV" for c in x) / len(x)),
}
rows = sorted(((k, spearmanr(v, ok.log2_expressed).statistic) for k, v in feats.items()),
              key=lambda kv: -abs(kv[1]))

fig, ax = plt.subplots(figsize=(10.5, 6.2))
names = [k for k, _ in rows]
r2 = [100 * r ** 2 for _, r in rows]
cols = [BLUE if "MonoYield" in n else SLATE for n in names]
y = np.arange(len(names))
ax.barh(y, r2, 0.64, color=cols)
for i, (v, (_, r)) in enumerate(zip(r2, rows)):
    ax.text(v + 0.04, i, f"{v:.1f}%", va="center", fontsize=14,
            color=BLUE if "MonoYield" in names[i] else DARK,
            fontweight="bold" if "MonoYield" in names[i] else "normal")
ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=13.5)
for lab, n in zip(ax.get_yticklabels(), names):
    if "MonoYield" in n:
        lab.set_color(BLUE); lab.set_fontweight("bold")
ax.invert_yaxis()
ax.set_xlim(0, 2.0)
ax.set_xlabel("rank variance in measured expression explained  "
              "(Spearman $\\rho^2$)", labelpad=10)

ax.set_title("A 6-billion-parameter yield model explains no more of display-sort\n"
             "expression than counting a protein's charged residues",
             fontsize=16, pad=16, color=DARK)
fig.text(0.5, 0.035,
         "5,114 de novo minibinders, all 85 aa — Chazin-Gray et al., Nat Commun 2026 · "
         "scored with ESMC-MonoYield (Setsu et al.)",
         ha="center", fontsize=11.5, color=DARK)
fig.text(0.5, 0.005,
         "Display-sort enrichment, not the soluble mg/L MonoYield predicts: the proxy is what fails here, not the model.",
         ha="center", fontsize=11.5, color=DARK)
fig.tight_layout(rect=[0, 0.075, 1, 1])
fig.savefig("figures/fimh_label_signal.png", dpi=220)
print("wrote figures/fimh_label_signal.png")
