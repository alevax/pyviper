import numpy as np
import pandas as pd
from scipy.stats import norm
from itertools import combinations
from typing import Dict, Optional

def ShadowRegulon_py(
    ss: pd.Series,
    nes: pd.Series,
    regul: Dict[str, Dict[str, pd.Series]],
    regulators: float = 0.05,
    shadow: float = 0.05,
    targets: int = 10,
    penalty: float = 2.0,
    method: str = "absolute"
) -> Optional[Dict[str, Dict[str, pd.Series]]]:
    # ---- method check (R: match.arg) ----
    method = str(method).lower()
    if method not in ("absolute", "adaptive"):
        raise ValueError("method must be 'absolute' or 'adaptive'")

    # ---- R-like coercions/names ----
    if not isinstance(ss, pd.Series):
        ss = pd.Series(ss)
    if not isinstance(nes, pd.Series):
        nes = pd.Series(nes)
    ss = ss.copy(); nes = nes.copy()
    ss.index = ss.index.astype(str)
    nes.index = nes.index.astype(str)

    # ---- two-sided p from NES under N(0,1) ----
    pval_all = pd.Series(2.0 * norm.sf(np.abs(nes.astype("float64"))), index=nes.index)

    # ---- choose candidate TFs: cutoff or top-N ----
    if regulators < 1:
        tfs = list(pval_all.index[pval_all.values < regulators])
    else:
        k = min(int(regulators), len(pval_all))
        order_idx = np.argsort(pval_all.values)[:k]
        tfs = list(pval_all.index[order_idx])

    # ---- drop composite names and require >=2 TFs ----
    tfs = [tf for tf in tfs if "--" not in tf]
    if len(tfs) < 2:
        return None

    # ---- unique(tfs) (order-preserving) ----
    unique_tfs = list(dict.fromkeys(tfs))

    # ---- per-tf1 list of per-tf2 one-sided p-values ----
    tmp: Dict[str, pd.Series] = {}

    for tf1 in unique_tfs:
        missing_tf1 = tf1 not in regul or 'tfmode' not in regul[tf1] or 'likelihood' not in regul[tf1]
        if not missing_tf1:
            tf1_tfmode = regul[tf1]['tfmode'].copy()
            tf1_like   = regul[tf1]['likelihood'].copy()
            tf1_tfmode.index = tf1_tfmode.index.astype(str)
            tf1_like.index   = tf1_like.index.astype(str)
        else:
            tf1_tfmode = pd.Series(dtype="float64")
            tf1_like   = pd.Series(dtype="float64")

        tf2_list = [x for x in tfs if x != tf1]

        # build overlap slices as a list (to allow duplicate tf2 labels, like R list)
        reg_overlap_list = []
        for tf2 in tf2_list:
            if missing_tf1 or tf2 not in regul or 'tfmode' not in regul[tf2]:
                reg_overlap_list.append((tf2, {
                    'tfmode': pd.Series(dtype="float64"),
                    'likelihood': pd.Series(dtype="float64")
                }))
                continue
            tf2_tfmode = regul[tf2]['tfmode'].copy()
            tf2_tfmode.index = tf2_tfmode.index.astype(str)
            pos_mask = tf1_tfmode.index.isin(tf2_tfmode.index)
            reg_overlap_list.append((
                tf2,
                {
                    'tfmode': tf1_tfmode.loc[pos_mask],
                    'likelihood': tf1_like.loc[pos_mask]
                }
            ))

        # s1/s2 transforms on genes regulated by tf1
        tf1_genes = tf1_tfmode.index
        ss_pos_mask = ss.index.isin(tf1_genes)
        ss_slice = ss.loc[ss_pos_mask].astype(float)

        if ss_slice.empty:
            tmp[tf1] = pd.Series([np.nan]*len(tf2_list), index=tf2_list, dtype="float64")
            continue

        ranks = ss_slice.rank(method="average")
        s2 = ranks / (len(ss_slice) + 1.0) * 2.0 - 1.0
        s1 = np.abs(s2) * 2.0 - 1.0
        if len(s1) > 0:
            s1 = s1 + (1.0 - float(s1.max())) / 2.0
        s1 = pd.Series(norm.ppf(s1 / 2.0 + 0.5), index=s1.index)

        sign_tf1 = float(np.sign(nes.get(tf1, 0.0))) or 1.0
        s2 = pd.Series(norm.ppf(s2 / 2.0 + 0.5), index=s2.index) * sign_tf1

        # per-pair stat → one-sided p (upper tail)
        per_pair_pvals = []
        for tf2, x in reg_overlap_list:
            x_tfmode = x['tfmode'].astype(float)
            x_like   = x['likelihood'].astype(float)

            if len(x_tfmode) < int(targets):
                per_pair_pvals.append(np.nan)
                continue

            s1_aligned = s1.reindex(x_tfmode.index)
            s2_aligned = s2.reindex(x_tfmode.index)

            sum1 = (x_tfmode * x_like * s2_aligned).sum(skipna=False)
            ss_sign = np.sign(sum1)
            if pd.isna(ss_sign):
                per_pair_pvals.append(np.nan)
                continue
            if ss_sign == 0:
                ss_sign = 1.0

            sum2 = ((1.0 - np.abs(x_tfmode)) * x_like * s1_aligned).sum(skipna=False)
            ww = x_like / x_like.max()
            denom = x_like.sum(skipna=False)
            core = (abs(sum1) + (sum2 if (not pd.isna(sum2) and sum2 > 0) else 0.0)) / denom
            stat = core * ss_sign * np.sqrt((ww**2).sum(skipna=False))

            per_pair_pvals.append(norm.sf(stat) if pd.notna(stat) else np.nan)

        tmp[tf1] = pd.Series(per_pair_pvals, index=tf2_list, dtype="float64")

    # ---- Flatten into a single named p-value vector ("A x B") ----
    flat_vals, flat_names = [], []
    for tf1 in unique_tfs:
        ser = tmp.get(tf1, pd.Series(dtype="float64"))
        tf2_names = list(ser.index)
        flat_vals.extend(ser.values)
        flat_names.extend([f"{tf1} x {tf2}" for tf2 in tf2_names])
    pval_vec = pd.Series(flat_vals, index=pd.Index(flat_names, dtype="object"), dtype="float64")

    # ---- Drop NAs ----
    pval_vec = pval_vec.dropna()

    # ---- Enumerate TF pairs and keep those computed ----
    regind = pd.DataFrame(list(combinations(tfs, 2)), columns=["A", "B"])
    pair_names = regind.apply(lambda r: f"{r['A']} x {r['B']}", axis=1)
    regind = regind[pair_names.isin(pval_vec.index)].reset_index(drop=True)

    # ---- Directional p-value matrix: p(A→B), p(B→A) ----
    ab_names = regind.apply(lambda r: f"{r['A']} x {r['B']}", axis=1)
    ba_names = regind.apply(lambda r: f"{r['B']} x {r['A']}", axis=1)
    pval_mat = pd.DataFrame({
        "A_to_B": pval_vec.reindex(ab_names).values,
        "B_to_A": pval_vec.reindex(ba_names).values
    }, index=regind.index)

    # ---- Count tests per TF (R: table(as.vector(regind))) ----
    tests = pd.Series(regind[["A", "B"]].values.ravel()).value_counts()

    # ==================================================================
    # R 'switch(method, absolute=..., adaptive=...)' — apply penalties
    # ==================================================================
    # Work on a copy so we don't mutate the input regul
    regul_corr: Dict[str, Dict[str, pd.Series]] = {
        tf: {
            "tfmode": (regul[tf]["tfmode"].copy() if "tfmode" in regul[tf] else pd.Series(dtype="float64")),
            "likelihood": (regul[tf]["likelihood"].copy() if "likelihood" in regul[tf] else pd.Series(dtype="float64")),
        }
        for tf in regul
    }
    for tf in regul_corr:
        regul_corr[tf]["tfmode"].index = regul_corr[tf]["tfmode"].index.astype(str)
        regul_corr[tf]["likelihood"].index = regul_corr[tf]["likelihood"].index.astype(str)

    changed: set = set()

    if method == "absolute":
        # Keep pairs where only one direction is significant (<shadow)
        mask_ab = (pval_mat["A_to_B"] < shadow) & (pval_mat["B_to_A"] > shadow)
        mask_ba = (pval_mat["A_to_B"] > shadow) & (pval_mat["B_to_A"] < shadow)
        tmp_abs = pd.concat(
            [
                regind.loc[mask_ab, ["A", "B"]].copy(),
                regind.loc[mask_ba, ["B", "A"]].rename(columns={"B": "A", "A": "B"})
            ],
            ignore_index=True
        )
        if tmp_abs.empty:
            return None

        # Penalize likelihood for the dominated TF (on overlapping targets)
        for _, row in tmp_abs.iterrows():
            dom, oth = row["A"], row["B"]
            if dom not in regul_corr or oth not in regul_corr:
                continue
            dom_tfmode = regul_corr[dom]["tfmode"]
            oth_tfmode = regul_corr[oth]["tfmode"]
            pos_mask = dom_tfmode.index.isin(oth_tfmode.index)

            ll = regul_corr[dom]["likelihood"].copy()
            expo = 1.0 / float(tests.get(dom, 1))  # 1 / tests[dom]
            ll.loc[pos_mask] = ll.loc[pos_mask] / (penalty ** expo)
            regul_corr[dom]["likelihood"] = ll
            changed.add(dom)

    else:  # adaptive
        # Directional evidence = log10(p(B→A)) - log10(p(A→B))
        pval1 = np.log10(pval_mat["B_to_A"].astype(float)) - np.log10(pval_mat["A_to_B"].astype(float))

        # Orient so column 1 is the dominated TF
        pos = pval1 > 0
        neg = pval1 < 0
        tmp_adp = pd.concat(
            [
                regind.loc[pos, ["A", "B"]].copy(),
                regind.loc[neg, ["B", "A"]].rename(columns={"B": "A", "A": "B"})
            ],
            ignore_index=True
        )
        # Magnitude always positive
        pval1_mag = pd.concat([pval1.loc[pos], -pval1.loc[neg]], ignore_index=True)

        if tmp_adp.empty:
            return None

        # Penalize likelihood proportionally to directional evidence
        for i, row in tmp_adp.iterrows():
            dom, oth = row["A"], row["B"]
            if dom not in regul_corr or oth not in regul_corr:
                continue
            w = float(pval1_mag.iloc[i])  # >= 0

            dom_tfmode = regul_corr[dom]["tfmode"]
            oth_tfmode = regul_corr[oth]["tfmode"]
            pos_mask = dom_tfmode.index.isin(oth_tfmode.index)

            ll = regul_corr[dom]["likelihood"].copy()
            expo = penalty / float(tests.get(dom, 1))  # penalty / tests[dom]
            ll.loc[pos_mask] = ll.loc[pos_mask] / ((1.0 + w) ** expo)
            regul_corr[dom]["likelihood"] = ll
            changed.add(dom)

    # ---- Return only those TFs that were actually modified ----
    if not changed:
        return None
    penalized_subset = {tf: regul_corr[tf] for tf in regul if tf in changed}
    return penalized_subset
