# Display-sort enrichment doesn't carry enough signal to validate a yield predictor

![](figures/fimh_label_signal.png)

5,114 de novo minibinders, every one 85 aa, each with a per-design expression readout —
including the ones that failed outright. Scored with
[ESMC-MonoYield](https://github.com/tsubo-lab/monoyield-predictor) using the authors' own
weights.

It explains **0.8%** of the rank variance. Counting charged residues gets **1.2%**, and those
two intervals overlap — nothing here is distinguishable from anything else, and nothing
reaches 2%. The limit is the label, not the model.

(The correlation is slightly negative, ρ = −0.09. At n = 5,021 that is significant,
p ≈ 1e-10, and still 0.8% of the variance. Across a different expression system it is not
evidence of anything directional, so it is a curiosity rather than a result.)

## The data

The expression readout was already public but not presented as one. Supplementary Data 5 of
[Chazin-Gray et al., *Nat Commun* (2026)](https://doi.org/10.1038/s41467-026-78123-6) is an NGS
pool table whose `sort` column includes an `expression` sort alongside `naive`, so every design
in the library has an expression measurement.

`data/fimh_display_expression.csv` — 5,114 rows. `log2_expressed` = log2(expression / naive),
spanning −7.2 to +4.8. `expression_flag` is `ok` (5,021), `zero_expression` (52, real naive
reads and zero expression reads) or `no_naive_reads` (41). Fixed length means length cannot
confound a sequence model.

`data/fimh_monoyield_predictions.csv` — ESMC-MonoYield output, mean of the five folds. It
is joined to the table above by row position, which the plotting script asserts against the
scored sequence rather than trusting.

## Caveats

**This is not the quantity MonoYield predicts, and not even the same expression system.**
ESMC-MonoYield is trained on monomeric yield in *E. coli*, measured in mg/L by a pooled
peptide-barcode + SEC + mass-spec assay under uniform conditions
([Setsu et al. 2026](https://doi.org/10.64898/2026.11.02.763507)). The readout here is NGS
enrichment from **cDNA display** (FimH, in vitro translation) and **yeast display** (Abp) —
a different organism, a different construct and a different quantity, usable rank-only
within this library. The table does not say which platform the 3,857 `Adhesin` rows used.

So this is a test of the label, not a test of the model. The authors' own prospective test is
the fair one, and it worked: **6 of 11** model-selected BindCraft designs exceeded 10 mg/L
monomeric yield, against **2 of 10** randomly selected. Nothing here contradicts that.

## Reproducing

```
python scripts/make_fimh_display_expression.py   # re-downloads the source table
python scripts/plot_fimh_label_signal.py         # rebuilds the figure
```

Data from the papers above, reused under their open-access terms. Code here: MIT.
