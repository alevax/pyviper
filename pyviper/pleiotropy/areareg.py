import numpy as np
import pandas as pd
from scipy.special import ndtri


def areareg(ss_i, sreg, minsize=5, mode="global"):
    ss_i = pd.to_numeric(pd.Series(ss_i, dtype=float), errors="coerce")
    ss_i.index = ss_i.index.astype(str)

    # precompute global transforms once
    r_all  = ss_i.rank(method="average")
    t2_all = r_all / (len(r_all) + 1.0)                # (0,1)
    t1_all = (t2_all - 0.5).abs() * 2.0
    if len(t1_all) > 0:
        t1_all = t1_all + (1.0 - float(t1_all.max())) / 2.0
    z2_all = pd.Series(ndtri(t2_all), index=ss_i.index)
    z1_all = pd.Series(ndtri(t1_all), index=ss_i.index)

    out = {}
    for tf, obj in sreg.items():
        m = pd.to_numeric(pd.Series(obj["tfmode"].values, index=obj["tfmode"].index.astype(str)), errors="coerce")
        l = pd.to_numeric(pd.Series(obj["likelihood"].values, index=obj["likelihood"].index.astype(str)), errors="coerce")

        keep = m.index.intersection(ss_i.index)
        if keep.size < minsize:
            continue

        m = m.loc[keep].values
        l = l.loc[keep].values

        if mode == "global":
            z2 = z2_all.loc[keep].values
            z1 = z1_all.loc[keep].values
        elif mode == "subset":
            sub = ss_i.loc[keep]
            r = sub.rank(method="average")
            t2 = r / (len(r) + 1.0)
            t1 = (t2 - 0.5).abs() * 2.0
            if len(t1) > 0:
                t1 = t1 + (1.0 - float(t1.max())) / 2.0
            z2 = ndtri(t2)
            z1 = ndtri(t1)
        else:
            raise ValueError("mode must be 'global' or 'subset'")

        if not np.isfinite(l).any():
            continue
        lmax = np.nanmax(np.abs(l))
        if not np.isfinite(lmax) or lmax == 0:
            continue

        s1 = np.sum(m * l * z2)
        sign = 1.0 if (np.isfinite(s1) and s1 >= 0) else -1.0
        s2 = np.sum((1.0 - np.abs(m)) * l * z1)

        denom = np.sum(l) + 1e-12
        ww = l / lmax
        scale = np.sqrt(np.sum(ww ** 2))

        core = (abs(s1) + (s2 if (np.isfinite(s2) and s2 > 0) else 0.0)) / denom
        out[str(tf)] = float(core * sign * scale)

    return pd.Series(out, dtype=float)
