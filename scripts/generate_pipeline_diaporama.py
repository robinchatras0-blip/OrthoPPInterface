import os
import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import imageio

# 1. Headless PyMOL setup
os.environ["PYOPENGL_PLATFORM"] = "osmesa"
import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()
from pymol import cmd

os.makedirs("wiki/slides", exist_ok=True)
os.makedirs("wiki/figures", exist_ok=True)

# Colors
DARK_BG = [0.043, 0.059, 0.098] # #0b0f19
HEX_DARK_BG = (11, 15, 25)
HEX_CYAN = (6, 182, 212)
HEX_AMBER = (245, 158, 11)
HEX_EMERALD = (16, 185, 129)
HEX_ROSE = (244, 63, 94)
HEX_BLUE = (56, 189, 248)
HEX_WHITE = (249, 250, 251)
HEX_MUTED = (156, 163, 175)

# Fonts
try:
    font_title = ImageFont.truetype("C:/Windows/Fonts/segoeuib.ttf", 28)
    font_subtitle = ImageFont.truetype("C:/Windows/Fonts/segoeui.ttf", 16)
    font_badge = ImageFont.truetype("C:/Windows/Fonts/consola.ttf", 14)
    font_hud_title = ImageFont.truetype("C:/Windows/Fonts/segoeuib.ttf", 15)
    font_hud_text = ImageFont.truetype("C:/Windows/Fonts/segoeui.ttf", 13)
    font_hud_mono = ImageFont.truetype("C:/Windows/Fonts/consola.ttf", 12)
except Exception:
    font_title = ImageFont.load_default()
    font_subtitle = font_title
    font_badge = font_title
    font_hud_title = font_title
    font_hud_text = font_title
    font_hud_mono = font_title

def setup_pymol_scene():
    cmd.reinitialize()
    cmd.set("bg_rgb", DARK_BG)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("stick_radius", 0.22)
    cmd.set("stick_ball", 1)
    cmd.set("stick_ball_ratio", 1.25)
    cmd.set("ray_shadows", 0)
    cmd.set("depth_cue", 1)
    cmd.set("antialias", 2)
    cmd.set("specular", 0.45)
    cmd.set("shininess", 35)

# Reference camera setup
ref_pdb = "data/inputs/complex_S1_S2.pdb"
setup_pymol_scene()
cmd.load(ref_pdb, "ref")
cmd.select("iface_box", "(chain A and resi 1-45) or (chain B and resi 1-50)")
cmd.orient("iface_box")
cmd.turn("x", -15)
cmd.turn("y", 30)
cmd.zoom("iface_box", 7)
saved_view = cmd.get_view()

# Mutations for cand_01 and cand_40
A_mut_resi = [3, 4, 12, 13, 14, 15, 16, 17, 18, 19, 29, 32, 33, 34, 37, 38, 39, 40, 41, 70, 71, 72, 185, 186, 201, 202, 203, 204, 205]
B_mut_resi = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19, 20, 29, 30, 32, 35, 37, 39, 40, 41, 42, 43, 44, 45, 47, 74, 75, 76, 77, 78, 79, 80, 81]

slides = [
    {
        "step_num": 1,
        "step_total": 6,
        "tag": "ÉTAT DE RÉFÉRENCE",
        "title": "Complexe Sauvage Endogène (A_WT · B_WT)",
        "subtitle": "Triade d'initiation ribosomique S1-S2 d'E. coli (Cible PDB : complexe 30S)",
        "script": "data/inputs/complex_S1_S2.pdb",
        "software": "Cristallographie PDB / BioPython DSSP",
        "params": "Interface contact < 6.0 Å (87 résidus sur A, 35 sur B). Enveloppe ARNr 16S à 4-6 Å.",
        "rationale": "Le recruteur S1 s'ancre obligatoirement sur S2 pour acheminer les ARNm vers le ribosome.",
        "mutations": "0 mutation (Sauvage natif) : S1_WT (224 res) · S2_WT (88 res)",
        "metrics": "Plafond sauvage empirique : iPTM = 0.5446 (interface masquée), ΔG = -48.2 REU",
        "pdb": "data/inputs/complex_S1_S2.pdb",
        "action": "wt"
    },
    {
        "step_num": 2,
        "step_total": 6,
        "tag": "MODULE 1 : CARTOGRAPHIE & TOPOLOGIE",
        "title": "Analyse d'Interface, Consolidation des Boucles & Origine Hotspots",
        "subtitle": "Prétraitement géométrique et définition des contraintes pour la diffusion",
        "script": "pipeline/src/01_analyze_interface.py",
        "software": "KDTree Distance Mapping / Algorithme SVD Kabsch",
        "params": "contact_cutoff: 6.0 Å, crop_distance: 15.0 Å, min_fixed_segment: 3, ppi_hotspots: 8",
        "rationale": "La fusion des îlots fixes courts (<3 res) en boucles continues empêche de casser le squelette par RFD3.",
        "mutations": "Délimitation de 3 boucles consolidées sur A (resi 3-22, 29-47, 185-205). Origine décalée de +10 Å.",
        "metrics": "87 résidus redesignables sur A, rognage de B à 15 Å (VRAM 24 Go GPU respectée)",
        "pdb": "data/inputs/complex_S1_S2.pdb",
        "action": "m1"
    },
    {
        "step_num": 3,
        "step_total": 6,
        "tag": "MODULE 2 : RUPTURE DE A'",
        "title": "Génération de la Rupture A' en État Apo",
        "subtitle": "Conception de novo de A' isolé pour abolir toute affinité envers B_WT",
        "script": "pipeline/src/02_generate_A_prime.py",
        "software": "RFdiffusion3 (apo-state) + LigandMPNN (rupture bias) + PyRosetta Ref2015",
        "params": "rupture_apo_design: true, mpnn_temp: 0.1, rupture_bias: inversion charges / répulsion",
        "rationale": "Diffuser en état apo évite que LigandMPNN n'optimise l'amarrage contre le partenaire sauvage B_WT.",
        "mutations": "29 mutations introduites sur A' : remodelage de poche alanine (G14A, T17A, R18A, F29A) et inversion (M3E, R4E)",
        "metrics": "Chaîne A' autonome générée et relaxée all-atom sous PyRosetta Ref2015",
        "pdb": "runs/run_overnight_production/05_final_eval/energy/A_prime_candidate_01_B_cand_40/Ap.BWT.pdb",
        "action": "m2"
    },
    {
        "step_num": 4,
        "step_total": 6,
        "tag": "MODULE 3 : FILTRAGE FAIL-FAST",
        "title": "Validation du Repliement & Rupture Sauvage (A' · B_WT)",
        "subtitle": "Portes éliminatoires rapides avant d'allouer les ressources GPU pour B'",
        "script": "pipeline/src/03_filter_A_prime.py",
        "software": "RoseTTAFold-3 All-Atom (RF3 monomère & co-amarrage)",
        "params": "pLDDT_monomer >= 80.0, RMSD <= 2.0 Å, iPTM(A'·B_WT) <= 0.35, max_motifs: 4",
        "rationale": "Élimine immédiatement les candidats mal repliés ou conservant un amarrage croisé parasite.",
        "mutations": "Validation de la lignée A'_candidate_01 par clustering de Hamming au plus éloigné",
        "metrics": "A' monomère : pLDDT = 84.2, RMSD = 1.15 Å. Rupture native : iPTM(A'·B_WT) = 0.295 (<= 0.35 validé)",
        "pdb": "runs/run_overnight_production/05_final_eval/energy/A_prime_candidate_01_B_cand_40/Ap.BWT.pdb",
        "action": "m3"
    },
    {
        "step_num": 5,
        "step_total": 6,
        "tag": "MODULE 4 : CO-ADAPTATION DE B'",
        "title": "Recalage SVD & Sauvetage Compensatoire de B'",
        "subtitle": "Conception de B' avec biais de spécificité différentielle contre A' fixé",
        "script": "pipeline/src/04_generate_B_prime.py",
        "software": "SVD Kabsch frame realignment + RFdiffusion3 multi-scaffold + LigandMPNN + PyRosetta",
        "params": "scaffold_max_flex_rmsd: 4.5 Å, framework_rmsd <= 1.5 Å, specificity_by_difference: true",
        "rationale": "Le biais favorise les contacts complémentaires aux mutations de A' et pénalise ceux avec A_WT.",
        "mutations": "37 mutations introduites sur B' : insertion du verrou aromatique (S12F, A47F) et nouveau pont salin (D39R)",
        "metrics": "Recalage SVD rigide (~70 Å corrigé), dérive de boucle bridée (< 4.5 Å), packing PyRosetta",
        "pdb": "runs/run_overnight_production/04_rescue_design/A_prime_candidate_01_B_cand_40.pdb",
        "action": "m4"
    },
    {
        "step_num": 6,
        "step_total": 6,
        "tag": "MODULE 5 : ÉVALUATION FINALE & VALIDATION",
        "title": "Triomphe de l'Orthogonalité : Sauvetage Cognat (A' · B')",
        "subtitle": "Amarrage stéréospécifique parfait et contre-validation indépendante AlphaFold2-Multimer v3",
        "script": "pipeline/src/05_eval_final.py & validate_af2_all19.py",
        "software": "RoseTTAFold-3 All-Atom (MSA masquées) + PyRosetta FastRelax + AF2-Multimer v3",
        "params": "unified_msa_scope: interface, dG_filter <= -30 REU, DockQ >= 0.23",
        "rationale": "L'interface mutée reconstitue un verrou moléculaire étanche incapable d'interagir avec l'hôte.",
        "mutations": "Complexe orthogonal complet : S1' (29 muts) · S2' (37 muts). Zéro interaction croisée avec l'hôte.",
        "metrics": "iPTM rescue = 0.745, DockQ = 0.4812, ΔG = -40.7 REU, ΔiPTM = +0.494 (AF2 RMSD = 1.79 Å vs 13.26 Å sauvage)",
        "pdb": "runs/run_overnight_production/05_final_eval/complex_rescue/A_prime_candidate_01_B_cand_40/A_prime_candidate_01_B_cand_40_unrelaxed_rank_001.pdb",
        "action": "m5"
    }
]

WIDTH = 1600
HEIGHT = 900
slide_paths = []

for idx, s in enumerate(slides):
    setup_pymol_scene()
    print(f"Rendering Step {s['step_num']}/{s['step_total']}: {s['title']}...")

    cmd.load(ref_pdb, "ref_frame")
    cmd.load(s["pdb"], "mol")

    # Rigidly align Chain A of mol to Chain A of ref_frame so the view stays strictly stable
    cmd.align("mol and chain A", "ref_frame and chain A")
    cmd.delete("ref_frame")

    cmd.set_view(saved_view)
    cmd.show("cartoon", "mol")
    cmd.hide("lines", "mol")

    action = s["action"]
    if action == "wt":
        cmd.color("cyan", "chain A")
        cmd.color("brightorange", "chain B")
        cmd.select("iface", "byres ((chain A within 6.0 of chain B) or (chain B within 6.0 of chain A))")
        cmd.show("sticks", "iface")
        cmd.color("cyan", "iface and chain A and elem C")
        cmd.color("brightorange", "iface and chain B and elem C")

    elif action == "m1":
        cmd.color("slate", "chain A")
        cmd.color("gray40", "chain B")
        cmd.set("cartoon_transparency", 0.5, "chain B")
        # Consolidated loops
        cmd.select("loops", "chain A and (resi 3-22 or resi 29-47 or resi 185-205)")
        cmd.color("teal", "loops")
        # Hotspots
        cmd.select("hotspots", "chain A and resi 12+14+15+17+18+19+32+38")
        cmd.show("sticks", "hotspots")
        cmd.color("yellow", "hotspots and elem C")

    elif action == "m2":
        cmd.hide("everything", "chain B")
        cmd.color("cyan", "chain A")
        mut_sel = "+".join(str(r) for r in A_mut_resi)
        cmd.select("a_muts", f"chain A and resi {mut_sel}")
        cmd.show("sticks", "a_muts")
        cmd.color("magenta", "a_muts and elem C")
        # Alanine cavity residues
        cmd.select("ala_cavity", "chain A and resi 14+17+18+29")
        cmd.color("hotpink", "ala_cavity and elem C")

    elif action == "m3":
        cmd.color("cyan", "chain A")
        cmd.color("brightorange", "chain B")
        mut_sel = "+".join(str(r) for r in A_mut_resi)
        cmd.select("a_muts", f"chain A and resi {mut_sel}")
        cmd.show("sticks", "a_muts")
        cmd.color("magenta", "a_muts and elem C")
        cmd.select("b_wt_iface", "chain B within 5.0 of a_muts")
        cmd.show("sticks", "b_wt_iface")
        cmd.color("brightorange", "b_wt_iface and elem C")

    elif action == "m4":
        cmd.color("cyan", "chain A")
        cmd.color("warmpink", "chain B")
        b_mut_sel = "+".join(str(r) for r in B_mut_resi)
        cmd.select("b_muts", f"chain B and resi {b_mut_sel}")
        cmd.show("sticks", "b_muts")
        cmd.color("tv_yellow", "b_muts and elem C")
        # Aromatic latch
        cmd.select("latch", "chain B and resi 12+47")
        cmd.show("spheres", "latch")
        cmd.set("sphere_scale", 0.45, "latch")
        cmd.color("chartreuse", "latch and elem C")

    elif action == "m5":
        cmd.color("cyan", "chain A")
        cmd.color("tv_yellow", "chain B")
        cmd.select("iface_rescue", "byres ((chain A within 5.5 of chain B) or (chain B within 5.5 of chain A))")
        cmd.show("sticks", "iface_rescue")
        cmd.color("cyan", "iface_rescue and chain A and elem C")
        cmd.color("tv_yellow", "iface_rescue and chain B and elem C")
        # Highlight locking latch
        cmd.select("lock_B", "chain B and resi 12+39+47")
        cmd.show("spheres", "lock_B")
        cmd.set("sphere_scale", 0.4, "lock_B")
        cmd.color("chartreuse", "lock_B and elem C")

    # Hetero-atom colors
    cmd.color("red", "elem O and (visible and not elem C)")
    cmd.color("blue", "elem N and (visible and not elem C)")
    cmd.color("yellow", "elem S and (visible and not elem C)")

    # Render raw image
    raw_path = f"wiki/slides/slide_{idx+1}_raw.png"
    cmd.png(raw_path, width=WIDTH, height=HEIGHT, ray=0)

    # Load into PIL
    img = Image.open(raw_path).convert("RGBA")
    overlay = Image.new("RGBA", (WIDTH, HEIGHT), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    # --- TOP BANNER (Glassmorphic dark card) ---
    draw.rectangle([(0, 0), (WIDTH, 115)], fill=(11, 15, 25, 235))
    draw.line([(0, 115), (WIDTH, 115)], fill=HEX_CYAN + (255,), width=3)

    # Badge & Step counter
    draw.rectangle([(30, 16), (280, 42)], fill=(6, 182, 212, 40), outline=HEX_CYAN + (200,), width=1)
    draw.text((40, 20), f"ORThoPPINTERFACE  |  ÉTAPE {s['step_num']}/{s['step_total']}", fill=HEX_CYAN, font=font_badge)

    draw.rectangle([(300, 16), (620, 42)], fill=(245, 158, 11, 40), outline=HEX_AMBER + (200,), width=1)
    draw.text((310, 20), s["tag"], fill=HEX_AMBER, font=font_badge)

    # Step progress dots on top right
    dot_start_x = WIDTH - 260
    for d_idx in range(1, s["step_total"] + 1):
        cx = dot_start_x + (d_idx - 1) * 36
        cy = 30
        if d_idx == s["step_num"]:
            draw.ellipse([(cx - 10, cy - 10), (cx + 10, cy + 10)], fill=HEX_CYAN, outline=HEX_WHITE, width=2)
            draw.text((cx - 4, cy - 8), str(d_idx), fill=HEX_DARK_BG, font=font_badge)
        elif d_idx < s["step_num"]:
            draw.ellipse([(cx - 7, cy - 7), (cx + 7, cy + 7)], fill=HEX_EMERALD)
            draw.text((cx - 3, cy - 7), str(d_idx), fill=HEX_DARK_BG, font=font_hud_mono)
        else:
            draw.ellipse([(cx - 6, cy - 6), (cx + 6, cy + 6)], fill=(55, 65, 81), outline=(100, 116, 139), width=1)

    # Title & Subtitle
    draw.text((30, 48), s["title"], fill=HEX_WHITE, font=font_title)
    draw.text((30, 84), s["subtitle"], fill=HEX_MUTED, font=font_subtitle)

    # --- BOTTOM HUD (Technical & Biophysical Data Card) ---
    hud_y = HEIGHT - 200
    draw.rectangle([(25, hud_y), (WIDTH - 25, HEIGHT - 20)], fill=(17, 24, 39, 240), outline=(55, 65, 81, 220), width=1)
    draw.line([(25, hud_y), (WIDTH - 25, hud_y)], fill=HEX_BLUE + (255,), width=2)

    col1_x = 45
    col2_x = 550
    col3_x = 1080

    # Column 1: Moteurs & Scripts
    draw.text((col1_x, hud_y + 15), "⚙️ LOGICIELS & MOTEURS", fill=HEX_AMBER, font=font_hud_title)
    draw.text((col1_x, hud_y + 40), f"Moteur : {s['software']}", fill=HEX_WHITE, font=font_hud_text)
    draw.text((col1_x, hud_y + 65), f"Script  : {s['script']}", fill=HEX_CYAN, font=font_hud_mono)
    draw.text((col1_x, hud_y + 92), "Raisonnement Biophysique :", fill=HEX_MUTED, font=font_hud_text)
    draw.text((col1_x, hud_y + 115), s["rationale"], fill=(226, 232, 240), font=font_hud_text)

    # Column 2: Paramètres
    draw.text((col2_x, hud_y + 15), "🔧 PARAMÈTRES PIPELINE", fill=HEX_EMERALD, font=font_hud_title)
    draw.text((col2_x, hud_y + 40), s["params"], fill=(203, 213, 225), font=font_hud_text)
    draw.text((col2_x, hud_y + 92), "Structure 3D Chargée :", fill=HEX_MUTED, font=font_hud_text)
    draw.text((col2_x, hud_y + 115), s["pdb"], fill=HEX_CYAN, font=font_hud_mono)

    # Column 3: Mutations & Métriques
    draw.text((col3_x, hud_y + 15), "📊 MUTATIONS & MÉTRIQUES RÉELLES", fill=HEX_BLUE, font=font_hud_title)
    draw.text((col3_x, hud_y + 40), s["mutations"], fill=(241, 245, 249), font=font_hud_text)
    draw.text((col3_x, hud_y + 92), "Résultats Empiriques Mesurés :", fill=HEX_MUTED, font=font_hud_text)
    draw.text((col3_x, hud_y + 115), s["metrics"], fill=(52, 211, 153), font=font_hud_text)

    # Merge overlay with rendered image
    final_img = Image.alpha_composite(img, overlay).convert("RGB")
    final_path = f"wiki/slides/slide_{idx+1}.png"
    final_img.save(final_path, quality=95)
    slide_paths.append(final_path)
    print(f"  -> Saved {final_path}")

cmd.quit()

# 2. Build Animated GIF (1.8s per slide)
print("\nCreating animated GIF...")
gif_images = [Image.open(p) for p in slide_paths]
gif_path = "wiki/figures/pipeline_diaporama_animation.gif"
# Resize for smooth web preview (1280x720)
gif_frames = [im.resize((1280, 720), Image.Resampling.LANCZOS) for im in gif_images]
gif_frames[0].save(
    gif_path,
    save_all=True,
    append_images=gif_frames[1:],
    duration=2200, # 2.2s per slide
    loop=0,
    optimize=True
)
print(f"Saved animated GIF to {gif_path} ({os.path.getsize(gif_path)} bytes)")

# 3. Build MP4 Video (4.0s per slide = 24 seconds total at 30 fps with smooth cross-dissolve)
print("\nCreating MP4 video with smooth transitions...")
video_path = "wiki/figures/pipeline_diaporama_video.mp4"
fps = 30
hold_seconds = 3.5
fade_seconds = 0.5
hold_frames = int(hold_seconds * fps)
fade_frames = int(fade_seconds * fps)

writer = imageio.get_writer(video_path, fps=fps, codec="libx264", quality=8)

for i in range(len(slide_paths)):
    curr_frame = np.array(Image.open(slide_paths[i]).convert("RGB"))
    next_frame = np.array(Image.open(slide_paths[(i + 1) % len(slide_paths)]).convert("RGB"))

    # Hold slide
    for _ in range(hold_frames):
        writer.append_data(curr_frame)

    # Cross-dissolve fade to next slide
    for f in range(fade_frames):
        alpha = (f + 1) / float(fade_frames)
        blended = (curr_frame * (1.0 - alpha) + next_frame * alpha).astype(np.uint8)
        writer.append_data(blended)

writer.close()
print(f"Saved MP4 video to {video_path} ({os.path.getsize(video_path)} bytes)")
print("ALL ASSETS CREATED SUCCESSFULLY!")
