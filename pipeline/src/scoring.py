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


IPSAE_PAE_CUTOFF = 10.0                                                # literature default (Dunbrack et al. 2025)
_IPSAE_EMPTY = {"ipsae": NAN, "ipsae_d0chn": NAN, "ipsae_d0dom": NAN, "lis": NAN}


def ipsae_score(pae, chain_ids, pae_cutoff=IPSAE_PAE_CUTOFF):
    """Interface-only pTM-like score computed straight from a PAE matrix, no extra folding needed
    (Dunbrack et al. 2025, https://www.biorxiv.org/content/10.1101/2025.02.10.637595, code at
    https://github.com/DunbrackLab/IPSAE). ipTM averages the TM-score-like transform of the PAE over the WHOLE
    chain pair, so a design can look "confident" from a well-folded, non-interacting core alone; ipSAE only
    averages over residue pairs that are both cross-chain AND low-PAE, so it can't be inflated that way -
    directly relevant to the MSA-conservation bias we measured in ipTM (README: native ipTM 0.865 full MSA vs
    0.559 with the interface masked). Faithful port of the d0res/d0chn/d0dom/LIS formulas, restricted to the
    two-chain case this pipeline always scores (returns NaN for anything else, e.g. a monomer fold).
    `chain_ids`: per-token chain label, same order/length as the rows/cols of `pae` (e.g. RF3's
    confidences.json `token_chain_ids`).
    """
    chains = np.asarray(chain_ids)
    pae = np.asarray(pae, dtype=float)
    labels = sorted(set(chains.tolist()))
    if len(labels) != 2 or pae.ndim != 2 or pae.shape[0] != len(chains) or pae.shape[1] != len(chains):
        return dict(_IPSAE_EMPTY)

    def d0(n):
        n = float(n)
        return max(1.0, 1.24 * (n - 15) ** (1.0 / 3.0) - 1.8) if n > 27 else 1.0

    def d0_array(n):
        return np.maximum(1.0, 1.24 * (np.maximum(26.0, n) - 15) ** (1.0 / 3.0) - 1.8)

    def ptm(x, d):
        return 1.0 / (1.0 + (x / d) ** 2.0)

    def direction(chain1, chain2):
        mask1, mask2 = chains == chain1, chains == chain2
        rows = np.where(mask1)[0]
        valid = np.outer(mask1, mask2) & (pae < pae_cutoff)

        ptm_d0chn = ptm(pae, d0(mask1.sum() + mask2.sum()))
        byres_d0chn = np.array([ptm_d0chn[i, valid[i]].mean() if valid[i].any() else 0.0 for i in rows])

        residues_1 = {i for i in rows if valid[i].any()}
        residues_2 = {j for i in residues_1 for j in np.where(valid[i])[0]}
        ptm_d0dom = ptm(pae, d0(len(residues_1) + len(residues_2)))
        byres_d0dom = np.array([ptm_d0dom[i, valid[i]].mean() if valid[i].any() else 0.0 for i in rows])

        d0res_byres = d0_array(valid.sum(axis=1))
        byres_d0res = np.array([ptm(pae[i, valid[i]], d0res_byres[i]).mean() if valid[i].any() else 0.0 for i in rows])

        cross_pae = pae[np.outer(mask1, mask2)]
        cross_pae = cross_pae[cross_pae < 12]
        lis = float(((12 - cross_pae) / 12).mean()) if cross_pae.size else 0.0

        return (float(byres_d0res.max()) if byres_d0res.size else 0.0,
                float(byres_d0chn.max()) if byres_d0chn.size else 0.0,
                float(byres_d0dom.max()) if byres_d0dom.size else 0.0, lis)

    fwd, rev = direction(labels[0], labels[1]), direction(labels[1], labels[0])
    return {"ipsae": max(fwd[0], rev[0]), "ipsae_d0chn": max(fwd[1], rev[1]),
            "ipsae_d0dom": max(fwd[2], rev[2]), "lis": (fwd[3] + rev[3]) / 2.0}


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
           "spearman_ipsae_vs_iptm_margins": rho("f_ipsae", "f_ortho_iptm"),
           "spearman_ipsae_margin_vs_binding": rho("f_ipsae", "f_energy"),
           "margin_sign_agreement": sign, "disagreements": []}
    if len(both) >= min_n:
        r_i, r_e = both["f_iptm_rel"].rank(), both["f_energy"].rank()
        far = (r_i - r_e).abs()
        rep["disagreements"] = [str(df.loc[i, "design_id"]) for i in far[far >= len(both) / 2].index]
    return {k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in rep.items()}
