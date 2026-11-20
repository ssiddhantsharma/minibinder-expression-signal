"""Per-design display-expression readouts for 5,114 de novo minibinders.

Source: FimH/Abp miniprotein inhibitor paper, Nat Commun 2026,
doi 10.1038/s41467-026-78123-6, Supplementary Data MOESM5 (the NGS pool table).
Its `sort` column carries an `expression` sort alongside `naive`, so every design
has an expression readout -- including the ones that failed outright.

NOT E. coli soluble yield: this is display expression, a proxy. The paper screened
FimH by cDNA display and Abp by yeast display, and the table does not say which
platform the `Adhesin` pools used.

    python scripts/make_fimh_display_expression.py
"""

import io
import urllib.request

import numpy as np
import pandas as pd

URL = ("https://static-content.springer.com/esm/"
       "art%3A10.1038%2Fs41467-026-78123-6/MediaObjects/"
       "41467_2026_78123_MOESM5_ESM.csv")
OUT = "data/fimh_display_expression.csv"


def main() -> None:
    raw = urllib.request.urlopen(URL, timeout=300).read()
    m = pd.read_csv(io.BytesIO(raw))
    n_in = m.sequence.nunique()

    piv = m.pivot_table(index="sequence", columns="sort", values="norm_count_frac", aggfunc="sum")
    meta = (m.sort_values("count", ascending=False).drop_duplicates("sequence")
            .set_index("sequence")[["description", "target", "pool_name"]])
    d = piv.join(meta, how="left").reset_index()

    d["binder_length"] = d.sequence.str.len()
    d["expressed_ratio"] = d["expression"] / d["naive"]
    # the edge cases are reported, never dropped: a design with zero expression reads
    # is the informative failure, not a row to discard
    d["expression_flag"] = np.where(d.naive.fillna(0) == 0, "no_naive_reads",
                            np.where(d["expression"].fillna(0) == 0, "zero_expression", "ok"))
    ok = d.expression_flag == "ok"
    d.loc[ok, "log2_expressed"] = np.log2(d.loc[ok, "expressed_ratio"])

    d["assay"] = "display_expression_sort_NGS"
    d["is_de_novo"] = True
    d["source"] = "BacterialAdhesins2025"
    d["source_url"] = "https://doi.org/10.1038/s41467-026-78123-6"

    assert len(d) == n_in, f"row count changed: {n_in} -> {len(d)}"
    assert d.sequence.duplicated().sum() == 0
    cols = ["sequence", "description", "target", "binder_length", "naive", "expression",
            "sort1", "sort2", "sort3", "expressed_ratio", "log2_expressed",
            "expression_flag", "assay", "is_de_novo", "source", "source_url", "pool_name"]
    d[cols].to_csv(OUT, index=False)
    print(f"wrote {OUT}  {len(d)} rows")
    print(" flags:", dict(d.expression_flag.value_counts()))


if __name__ == "__main__":
    main()
