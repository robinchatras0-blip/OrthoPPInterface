import os
import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

# 1. Headless PyMOL setup
os.environ["PYOPENGL_PLATFORM"] = "osmesa"
import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()
from pymol import cmd

os.makedirs("wiki/slides", exist_ok=True)
os.makedirs("wiki/figures", exist_ok=True)

# Common visual styles
DARK_BG = [0.043, 0.059, 0.098] # #0b0f19
CYAN_A = [0.06, 0.75, 0.85]     # Chain A / A'
ORANGE_B = [0.95, 0.60, 0.15]   # Chain B_WT
GOLD_BP = [0.95, 0.75, 0.10]    # Chain B' rescue
MAGENTA_MUT = [0.95, 0.25, 0.60]# A' mutations
LIME_MUT = [0.20, 0.85, 0.40]   # B' mutations

def setup_scene():
    cmd.reinitialize()
    cmd.set("bg_rgb", DARK_BG)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("stick_radius", 0.20)
    cmd.set("stick_ball", 1)
    cmd.set("stick_ball_ratio", 1.2)
    cmd.set("ray_shadows", 0)
    cmd.set("depth_cue", 1)
    cmd.set("antialias", 2)
    cmd.set("specular", 0.4)
    cmd.set("shininess", 30)

# Reference structure for shared camera orientation
ref_pdb = "data/inputs/complex_S1_S2.pdb"
setup_scene()
cmd.load(ref_pdb, "ref")
cmd.select("interface_core", "(chain A and resi 1-50) or (chain B and resi 1-50)")
cmd.orient("interface_core")
cmd.turn("x", -15)
cmd.turn("y", 25)
cmd.zoom("interface_core", 8)
saved_view = cmd.get_view()
print("Saved reference camera orientation successfully.")

# Exact mutation positions for cand_01 and cand_40
A_mut_resi = [3, 4, 12, 13, 14, 15, 16, 17, 18, 19, 29, 32, 33, 34, 37, 38, 39, 40, 41, 70, 71, 72, 185, 186, 201, 202, 203, 204, 205]
B_mut_resi = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19, 20, 29, 30, 32, 35, 37, 39, 40, 41, 42, 43, 44, 45, 47, 74, 75, 76, 77, 78, 79, 80, 81]

slides_config = [
    {
        "id": "slide1_native_complex",
        "title": "Étape 0 : Complexe Endogène Initial (A_WT · B_WT)",
        "module": "État de Référence Biologique",
        "script": "data/inputs/complex_S1_S2.pdb",
        "software": "BioPython PDBParser / PDB Crystallography",
        "params": "Structure S1-S2 30S cristallographique native (224 résidus A, 88 résidus B)",
        "metrics": "iPTM ceiling = 0.5446 (interface masquée), ΔG = -48.2 REU",
        "pdb": "data/inputs/complex_S1_S2.pdb",
        "action": "wt"
    },
    {
        "id": "slide2_module1_topology",
        "title": "Étape 1 : Cartographie Topologique & Consolidation des Boucles",
        "module": "Module 1 : Analyse de l'Interface",
        "script": "pipeline/src/01_analyze_interface.py",
        "software": "KDTree distance mapping, SVD coordinate geometry",
        "params": "contact_threshold: 6.0 Å, crop_distance: 15.0 Å, min_fixed_segment: 3, ppi_hotspots: 8",
        "metrics": "87 résidus interface sur A, 35 sur B. 3 boucles consolidées (20/19/8 res), origine +10 Å",
        "pdb": "data/inputs/complex_S1_S2.pdb",
        "action": "m1_topology"
    },
    {
        "id": "slide3_module2_a_prime_rupture",
        "title": "Étape 2 : Diffusion Apo & Rupture de Séquence sur A'",
        "module": "Module 2 : Génération de A'",
        "script": "pipeline/src/02_generate_A_prime.py",
        "software": "RFdiffusion3 (apo-state), LigandMPNN (rupture bias), PyRosetta FastRelax",
        "params": "rupture_apo_design: true, mpnn_temp: 0.1, rupture_bias: inversion charges/stérique",
        "metrics": "29 mutations introduites sur A' (remodelage alanine des poches: G14A, T17A, R18A, F29A)",
        "pdb": "runs/run_overnight_production/05_final_eval/energy/A_prime_candidate_01_B_cand_40/Ap.BWT.pdb",
        "action": "m2_apo"
    },
    {
        "id": "slide4_module3_fail_fast_rupture",
        "title": "Étape 3 : Portes Fail-Fast & Rupture Sauvage Confirmée (A' · B_WT)",
        "module": "Module 3 : Filtrage Fail-Fast",
        "script": "pipeline/src/03_filter_A_prime.py",
        "software": "RoseTTAFold-3 All-Atom (RF3 inference mode)",
        "params": "pLDDT_cutoff: >=80.0, RMSD_cutoff: <=2.0 Å, rupture_iPTM_cutoff: <=0.35, max_motifs: 4",
        "metrics": "pLDDT(A') = 84.2, RMSD = 1.15 Å. Rupture native: iPTM(A'·B_WT) = 0.295 (seuil <= 0.35 validé)",
        "pdb": "runs/run_overnight_production/05_final_eval/energy/A_prime_candidate_01_B_cand_40/Ap.BWT.pdb",
        "action": "m3_rupture"
    },
    {
        "id": "slide5_module4_b_prime_rescue",
        "title": "Étape 4 : Recalage SVD & Co-Adaptation Sauvetage sur B'",
        "module": "Module 4 : Génération de B'",
        "script": "pipeline/src/04_generate_B_prime.py",
        "software": "SVD Kabsch frame realignment, RFdiffusion3 multi-scaffold, LigandMPNN spécifique, PyRosetta",
        "params": "scaffold_max_flex_rmsd: 4.5 Å, framework_rmsd: <=1.5 Å, specificity_by_difference: true",
        "metrics": "37 mutations compensatoires sur B' (S12F, A47F verrou aromatique, D39R nouveau pont salin)",
        "pdb": "runs/run_overnight_production/04_rescue_design/A_prime_candidate_01_B_cand_40.pdb",
        "action": "m4_rescue"
    },
    {
        "id": "slide6_module5_cognate_champions",
        "title": "Étape 5 : Complexe Sauvé Cognat (A' · B') & Orthogonalité Parfaite",
        "module": "Module 5 : Évaluation Finale & Validation",
        "script": "pipeline/src/05_eval_final.py & validate_af2_all19.py",
        "software": "RoseTTAFold-3 All-Atom (masked MSA), PyRosetta Ref2015, DockQ CAPRI, AF2-Multimer v3",
        "params": "unified_msa_scope: interface, dG_filter: <= -30 REU, DockQ_gate: >= 0.23",
        "metrics": "iPTM rescue = 0.745, DockQ = 0.4812, dG = -40.7 REU, ΔiPTM = +0.494 (AF2 RMSD = 1.79 Å)",
        "pdb": "runs/run_overnight_production/05_final_eval/complex_rescue/A_prime_candidate_01_B_cand_40/A_prime_candidate_01_B_cand_40_unrelaxed_rank_001.pdb",
        "action": "m5_final"
    }
]

rendered_images = []

for idx, sc in enumerate(slides_config):
    setup_scene()
    print(f"Rendering slide {idx+1}/{len(slides_config)}: {sc['id']}...")
    cmd.load(sc["pdb"], "mol")
    cmd.set_view(saved_view)

    # General styling
    cmd.show("cartoon", "mol")
    cmd.hide("lines", "mol")

    action = sc["action"]
    if action == "wt":
        cmd.color("cyan", "chain A")
        cmd.color("brightorange", "chain B")
        cmd.select("iface", "byres ((chain A within 6.0 of chain B) or (chain B within 6.0 of chain A))")
        cmd.show("sticks", "iface")
        cmd.color("cyan", "iface and chain A and elem C")
        cmd.color("brightorange", "iface and chain B and elem C")

    elif action == "m1_topology":
        cmd.color("slate", "chain A")
        cmd.color("gray40", "chain B")
        # Highlight consolidated loop regions
        cmd.select("loops", "chain A and (resi 3-22 or resi 29-47 or resi 185-205)")
        cmd.color("teal", "loops")
        # Highlight hotspots in sticks
        cmd.select("hotspots", "chain A and resi 12+14+15+17+18+19+32+38")
        cmd.show("sticks", "hotspots")
        cmd.color("yellow", "hotspots and elem C")

    elif action == "m2_apo":
        cmd.hide("everything", "chain B")
        cmd.color("cyan", "chain A")
        # Highlight 29 mutated positions in magenta sticks
        mut_sel = "+".join(str(r) for r in A_mut_resi)
        cmd.select("a_muts", f"chain A and resi {mut_sel}")
        cmd.show("sticks", "a_muts")
        cmd.color("magenta", "a_muts and elem C")

    elif action == "m3_rupture":
        cmd.color("cyan", "chain A")
        cmd.color("brightorange", "chain B")
        mut_sel = "+".join(str(r) for r in A_mut_resi)
        cmd.select("a_muts", f"chain A and resi {mut_sel}")
        cmd.show("sticks", "a_muts")
        cmd.color("magenta", "a_muts and elem C")
        cmd.select("b_wt_iface", "chain B within 5.0 of a_muts")
        cmd.show("sticks", "b_wt_iface")
        cmd.color("brightorange", "b_wt_iface and elem C")

    elif action == "m4_rescue":
        cmd.color("cyan", "chain A")
        cmd.color("warmpink", "chain B")
        b_mut_sel = "+".join(str(r) for r in B_mut_resi)
        cmd.select("b_muts", f"chain B and resi {b_mut_sel}")
        cmd.show("sticks", "b_muts")
        cmd.color("tv_yellow", "b_muts and elem C")
        # Highlight key aromatic latch
        cmd.select("latch", "chain B and resi 12+47")
        cmd.show("spheres", "latch")
        cmd.set("sphere_scale", 0.45, "latch")

    elif action == "m5_final":
        cmd.color("cyan", "chain A")
        cmd.color("tv_yellow", "chain B")
        cmd.select("iface_rescue", "byres ((chain A within 5.5 of chain B) or (chain B within 5.5 of chain A))")
        cmd.show("sticks", "iface_rescue")
        cmd.color("cyan", "iface_rescue and chain A and elem C")
        cmd.color("tv_yellow", "iface_rescue and chain B and elem C")
        # Highlight specific locking residues
        cmd.select("lock_A", "chain A and resi 14+17+18")
        cmd.select("lock_B", "chain B and resi 12+39+47")
        cmd.show("spheres", "lock_B")
        cmd.set("sphere_scale", 0.4, "lock_B")

    # Standard element coloring for all sticks
    cmd.color("red", "elem O and (visible and not elem C)")
    cmd.color("blue", "elem N and (visible and not elem C)")
    cmd.color("yellow", "elem S and (visible and not elem C)")

    # Save raw render image
    raw_img_path = f"wiki/slides/{sc['id']}_raw.png"
    cmd.png(raw_img_path, width=1280, height=720, ray=0)
    
    # Post-process with PIL to add HUD / info overlay
    pil_img = Image.open(raw_img_path)
    draw = ImageDraw.Draw(pil_img)
    
    # Overlay banner at top
    draw.rectangle([(0, 0), (1280, 85)], fill=(11, 15, 25, 220))
    draw.line([(0, 85), (1280, 85)], fill=(6, 182, 212), width=2)
    
    # Overlay banner at bottom for parameters
    draw.rectangle([(0, 570), (1280, 720)], fill=(11, 15, 25, 230))
    draw.line([(0, 570), (1280, 570)], fill=(56, 189, 248), width=2)

    # Use default font
    # Title
    draw.text((25, 12), f"Étape {idx+1} / 6 : {sc['title']}", fill=(249, 250, 251))
    draw.text((25, 36), f"Module : {sc['module']}  |  Script : {sc['script']}", fill=(6, 182, 212))
    draw.text((25, 58), f"Champion Benchmark : A_prime_candidate_01 · B_cand_40 (S1-S2 Ribosomal)", fill=(156, 163, 175))

    # Bottom parameters
    draw.text((25, 580), "⚙️ LOGICIELS & MOTEURS :", fill=(245, 158, 11))
    draw.text((230, 580), sc["software"], fill=(226, 232, 240))

    draw.text((25, 605), "🔧 PARAMÈTRES PIPELINE :", fill=(16, 185, 129))
    draw.text((230, 605), sc["params"], fill=(203, 213, 225))

    draw.text((25, 630), "📊 RÉSULTATS & MÉTRIQUES :", fill=(56, 189, 248))
    draw.text((230, 630), sc["metrics"], fill=(241, 245, 249))

    final_img_path = f"wiki/slides/{sc['id']}.png"
    pil_img.save(final_img_path)
    rendered_images.append(final_img_path)
    print(f"Saved slide to {final_img_path}")

cmd.quit()
print(f"All {len(rendered_images)} slides generated successfully!")
