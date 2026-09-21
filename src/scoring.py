"""Orthogonality scores: RF3-based, energy-based and combined. Pure functions (no GPU / PyRosetta needed).

Both margins have the same shape "rescue - best wild-type cross-reaction" and are normalised so that the native complex
is 1, which makes them combinable:
    F_iptm   = (iPTM(A'.B') - max(iPTM(A'.B_WT), iPTM(A_WT.B'))) / iPTM_ceiling         (ceiling = WT/WT, same MSA masks)
    F_energy = b(A'.B') - max(b(A'.B_WT), b(A_WT.B')),   b = dG / dG_native              (dG < 0 = binding)
    F_ortho  = w * F_iptm + (1 - w) * F_energy
"""
import math

import numpy as np

NAN = float('nan')


def ok(*xs):
    return all(x is not None and x == x for x in xs)


def iptm_margin(iptm_rescue, iptm_rupture, iptm_negative):
    """Raw RF3 margin (the original F_ortho)."""
    return iptm_rescue - max(iptm_rupture, iptm_negative) if ok(iptm_rescue, iptm_rupture, iptm_negative) else NAN


def f_iptm_rel(iptm_rescue, iptm_rupture, iptm_negative, ceiling):
    m = iptm_margin(iptm_rescue, iptm_rupture, iptm_negative)
    return m / ceiling if ok(m, ceiling) and ceiling > 0 else NAN


def energy_fields(dG_rescue, dG_negative, dG_rupture, dG_native):
    """Binding fractions relative to the native complex, cross-pair energy gaps and the energy margin."""
    empty = {"b_rescue": NAN, "b_negative": NAN, "b_rupture": NAN, "f_energy": NAN, "gap_negative": NAN, "gap_rupture": NAN}
    if not ok(dG_native) or dG_native >= 0:
        return empty
    b = {k: (v / dG_native if ok(v) else NAN) for k, v in (("rescue", dG_rescue), ("negative", dG_negative), ("rupture", dG_rupture))}
    out = {"b_rescue": b["rescue"], "b_negative": b["negative"], "b_rupture": b["rupture"],
           "gap_negative": dG_rescue - dG_negative if ok(dG_rescue, dG_negative) else NAN,
           "gap_rupture": dG_rescue - dG_rupture if ok(dG_rescue, dG_rupture) else NAN}
    out["f_energy"] = b["rescue"] - max(b["negative"], b["rupture"]) if ok(*b.values()) else NAN
    return out


def combined_f_ortho(f_iptm, f_energy, weight_iptm, energy_enabled):
    """Weighted combination; iPTM only when the energy stage is disabled. NaN if a needed component is missing."""
    if not energy_enabled:
        return f_iptm
    return weight_iptm * f_iptm + (1.0 - weight_iptm) * f_energy if ok(f_iptm, f_energy) else NAN


def coherence_report(df, min_n=4):
    """How well RF3 and the energy agree across designs (Spearman on ranks; sign agreement of the two margins)."""
    def rho(a, b):
        sub = df[[a, b]].dropna() if a in df and b in df else df.iloc[0:0]
        # Spearman = Pearson on ranks (computed here: pandas needs scipy for method="spearman")
        if len(sub) < min_n:
            return NAN
        with np.errstate(invalid="ignore", divide="ignore"):                      # constant column -> NaN
            return float(sub[a].rank().corr(sub[b].rank()))

    both = df[["f_iptm_rel", "f_energy"]].dropna() if "f_iptm_rel" in df and "f_energy" in df else df.iloc[0:0]
    sign = float(((both["f_iptm_rel"] > 0) == (both["f_energy"] > 0)).mean()) if len(both) else NAN
    rep = {"n_designs": int(len(df)), "n_with_both": int(len(both)),
           "spearman_rescue_iptm_vs_binding": rho("iptm_rescue", "b_rescue"),
           "spearman_negative_iptm_vs_binding": rho("iptm_negative", "b_negative"),
           "spearman_rupture_iptm_vs_binding": rho("iptm_rupture", "b_rupture"),
           "spearman_margins": rho("f_iptm_rel", "f_energy"),
           "margin_sign_agreement": sign, "disagreements": []}
    if len(both) >= min_n:
        r_i, r_e = both["f_iptm_rel"].rank(), both["f_energy"].rank()
        far = (r_i - r_e).abs()
        rep["disagreements"] = [str(df.loc[i, "design_id"]) for i in far[far >= len(both) / 2].index]
    return {k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in rep.items()}
