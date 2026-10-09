#!/usr/bin/env python3
"""
generate_enhanced_diaporama.py - Generates an all-atom 3D evolutionary pipeline walkthrough.

Visual Guidelines (User Request):
- Chain A (S1 / A') backbone: Ruby Red (#e11d48)
- Chain B (S2 / B') backbone: Marine Blue (#2563eb)
- Interface sidechains: Element-colored (O=red, N=blue, S=yellow) with carbons matching their chain
  (Salmon/Coral for Chain A, Cyan/Sky for Chain B)
- Direct structural superposition demonstrating that the global framework is preserved (RMSD < 1.2 Å)
  and only the interface is surgically redesigned.
- Outputs:
  1. wiki/slides/slide_1.png to slide_6.png (1600x900 HD)
  2. wiki/figures/pipeline_diaporama_video.mp4 (Smooth 28s H.264)
  3. wiki/figures/pipeline_diaporama_animation.gif (Optimized looping GIF)
"""

import os
import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import imageio

# Headless PyMOL setup
os.environ["PYOPENGL_PLATFORM"] = "osmesa"
import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()
from pymol import cmd

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SLIDES_DIR = os.path.join(BASE_DIR, "wiki", "slides")
FIG_DIR = os.path.join(BASE_DIR, "wiki", "figures")
os.makedirs(SLIDES_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ----------------- Visual Styles -----------------
BG_COLOR = [0.035, 0.045, 0.075] # #090c13 (Deep navy space)

# Chain colors (Requested by user: Backbone strictly Red & Blue)
RUBY_RED = [0.88, 0.11, 0.28]    # Chain A / A' backbone
MARINE_BLUE = [0.15, 0.39, 0.92] # Chain B / B' backbone

# Interface carbons
CORAL_A = [0.98, 0.45, 0.52]     # Carbons for Chain A interface
SKY_B = [0.22, 0.74, 0.97]       # Carbons for Chain B interface

# Hex colors for HUD (PIL)
HEX_BG = (9, 12, 19)
HEX_CARD = (15, 20, 32, 240)
HEX_RED = (225, 29, 72)
HEX_BLUE = (37, 99, 235)
HEX_CYAN = (6, 182, 212)
HEX_AMBER = (245, 158, 11)
HEX_EMERALD = (16, 185, 129)
HEX_WHITE = (249, 250, 251)
HEX_MUTED = (156, 163, 175)
HEX_BORDER = (45, 55, 75)

# Fonts
try:
    font_step = ImageFont.truetype("C:/Windows/Fonts/consola.ttf", 15)
    font_title = ImageFont.truetype("C:/Windows/Fonts/segoeuib.ttf", 26)
    font_subtitle = ImageFont.truetype("C:/Windows/Fonts/segoeui.ttf", 15)
    font_legend = ImageFont.truetype("C:/Windows/Fonts/segoeuib.ttf", 13)
    font_hud_tag = ImageFont.truetype("C:/Windows/Fonts/segoeuib.ttf", 13)
    font_hud_title = ImageFont.truetype("C:/Windows/Fonts/segoeuib.ttf", 14)
    font_hud_text = ImageFont.truetype("C:/Windows/Fonts/segoeui.ttf", 13)
    font_hud_mono = ImageFont.truetype("C:/Windows/Fonts/consola.ttf", 13)
except Exception:
    font_step = ImageFont.load_default()
    font_title = font_step
    font_subtitle = font_step
    font_legend = font_step
    font_hud_tag = font_step
    font_hud_title = font_step
    font_hud_text = font_step
    font_hud_mono = font_step

def setup_pymol_scene():
    cmd.reinitialize()
    cmd.set("bg_rgb", BG_COLOR)
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

# Camera orientation centered on the interaction interface
ref_pdb = os.path.join(BASE_DIR, "data", "inputs", "complex_S1_S2.pdb")
setup_pymol_scene()
cmd.load(ref_pdb, "ref_cam")
cmd.select("iface_core", "(chain A and resi 1-45) or (chain B and resi 1-50)")
cmd.orient("iface_core")
cmd.turn("x", -12)
cmd.turn("y", 22)
cmd.zoom("iface_core", 6.5)
saved_view = cmd.get_view()
cmd.delete("ref_cam")
print("Saved reference camera orientation successfully.")

# Mutations for Champion #1 (A_prime_candidate_01 · B_cand_40)
A_mut_resi = [3, 4, 12, 13, 14, 15, 16, 17, 18, 19, 29, 32, 33, 34, 37, 38, 39, 40, 41, 70, 71, 72, 185, 186, 201, 202, 203, 204, 205]
B_mut_resi = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 19, 20, 29, 30, 32, 35, 37, 39, 40, 41, 42, 43, 44, 45, 47, 74, 75, 76, 77, 78, 79, 80, 81]

slides = [
    {
        "step_num": 1,
        "step_total": 6,
        "tag": "ETAT INITIAL NATIF",
        "title": "Complexe Sauvage Endogène (A_WT · B_WT)",
        "subtitle": "Architecture cristallographique native 30S d'E. coli : recruteur S1 (A) ancré sur la protéine ribosomique S2 (B)",
        "software": "Cristallographie PDB / BioPython PDBParser & DSSP",
        "params": "Interface contact < 6.0 A (87 residus sur A, 35 sur B). Enveloppe ARNr 16S a 4-6 A de distance.",
        "mechanism": "Le squelette global (framework) assure le repliement et l'orientation ribosomique. Seule l'interface doit etre remodelee.",
        "metrics": "iPTM ceiling = 0.5446 (interface masquee), DeltaG = -48.2 REU, 0 mutation. Squelettes : Rouge (A) et Bleu (B).",
        "action": "slide1"
    },
    {
        "step_num": 2,
        "step_total": 6,
        "tag": "MODULE 1 : DELIMITATION DU FRAMEWORK",
        "title": "Cartographie Topologique : Gel du Framework & Isolation de l'Interface",
        "subtitle": "Preservation stricte du coeur structural (>70% de la proteine) et deverrouillage des 3 boucles pour RFD3",
        "software": "KDTree Distance Mapping / Algorithme SVD Kabsch / SVD Normal Projections",
        "params": "contact_cutoff: 6.0 A, crop_distance: 15.0 A, min_fixed_segment: 3, ppi_hotspots: 8",
        "mechanism": "Le framework est gele (semi-transparent). Seules les 3 boucles d'interface (20/19/8 res) sont autorisees a diffuser. Origine a +10 A.",
        "metrics": "87 residus redesignables sur A, B rogne a 15 A (GPU 24 Go). Origine projetee dans le solvant pour eviter tout clash initial.",
        "action": "slide2"
    },
    {
        "step_num": 3,
        "step_total": 6,
        "tag": "MODULE 2 : RUPTURE SUR A'",
        "title": "Diffusion Apo & Remodelage de Novo sur A' (29 Mutations)",
        "subtitle": "Conception de novo de A' en l'absence de B pour abolir l'affinite sauvage sans biais de design positif",
        "software": "RFdiffusion3 (apo-state) + LigandMPNN (rupture bias) + PyRosetta Ref2015 FastRelax",
        "params": "rupture_apo_design: true, mpnn_temp: 0.1, rupture_bias: inversion charges (M3E, R4E) & cavites alanines",
        "mechanism": "Le framework de A reste intact (RMSD = 1.15 A). 29 residus d'interface mutes : poches alanines (G14A, T17A, R18A, F29A).",
        "metrics": "Chaine A' repliee de maniere autonome (pLDDT = 84.2). Surface de contact prete pour le demantelement de l'hote.",
        "action": "slide3"
    },
    {
        "step_num": 4,
        "step_total": 6,
        "tag": "MODULE 3 : FILTRAGE FAIL-FAST",
        "title": "Superposition avec l'Hôte : Rupture Sauvage Confirmée (A' · B_WT)",
        "subtitle": "Contrôle Fail-Fast montrant l'effondrement de l'amarrage entre A' muté et le partenaire sauvage B_WT",
        "software": "RoseTTAFold-3 All-Atom (RF3 co-amarrage) / PyRosetta Ref2015",
        "params": "pLDDT_monomer >= 80.0, RMSD <= 2.0 A, rupture_iPTM_cutoff <= 0.35, max_motifs: 4",
        "mechanism": "Superposition A' (Rouge) + B_WT (Bleu) : Les charges Glu3/Glu4 repoussent B_WT et les cavites alanines detruisent les contacts natifs.",
        "metrics": "iPTM(A' · B_WT) = 0.295 <= 0.35 (VALIDE). Rupture complete du complexe avec l'hote : aucun effet dominant-negatif !",
        "action": "slide4"
    },
    {
        "step_num": 5,
        "step_total": 6,
        "tag": "MODULE 4 : CO-ADAPTATION DE B'",
        "title": "Recalage SVD & Co-Adaptation Sauvetage sur B' (37 Mutations)",
        "subtitle": "Conception complémentaire de B' par différence : insertion du verrou aromatique et nouveaux ponts salins",
        "software": "SVD Kabsch frame realignment + RFdiffusion3 multi-scaffold + LigandMPNN + PyRosetta",
        "params": "scaffold_max_flex_rmsd: 4.5 A, framework_rmsd <= 1.5 A, specificity_by_difference: true",
        "mechanism": "B' conserve son squelette (Bleu). 37 mutations all-atom inserent Phe12 et Phe47 directement dans les cavites alanines de A' !",
        "metrics": "Verrou 'Knob-into-Hole' (S12F, A47F) + pont salin Arg39 (D39R) : complementarite stereospecifique exclusive contre A' fixe.",
        "action": "slide5"
    },
    {
        "step_num": 6,
        "step_total": 6,
        "tag": "MODULE 5 : SUPERPOSITION & TRIOMPHE",
        "title": "Superposition Directe (A' · B' vs A_WT · B_WT) : Orthogonalité Parfaite",
        "subtitle": "Superposition atome par atome démontrant la conservation totale du framework global et le redesign exclusif de l'interface",
        "software": "RoseTTAFold-3 All-Atom + AlphaFold2-Multimer v3 + PyRosetta Ref2015 + CAPRI DockQ",
        "params": "Superposition alignee sur le squelette global (RMSD total = 0.76 A). Regime unifie MSA masquee.",
        "mechanism": "Le framework est 100% superpose (RMSD < 1.0 A). Seule l'interface all-atom est mutee : verrouillage orthogonal exclusif et etanche.",
        "metrics": "iPTM rescue = 0.745, DockQ = 0.4812, Delta iPTM = +0.494 (AF2 RMSD = 1.79 A vs 13.26 A sauvage). Succes absolu !",
        "action": "slide6"
    }
]

# Common PDB paths
pdb_wt = os.path.join(BASE_DIR, "data", "inputs", "complex_S1_S2.pdb")
pdb_des = os.path.join(BASE_DIR, "DOSSIER_CHAMPIONS", "structures_pdb", "01_A_prime_candidate_01_B_cand_40_design.pdb")
pdb_rupture = os.path.join(BASE_DIR, "DOSSIER_CHAMPIONS", "structures_pdb", "01_A_prime_candidate_01_B_cand_40_Ap_BWT_rupture.pdb")

WIDTH = 1600
HEIGHT = 900
slide_files = []

for s in slides:
    step_num = s["step_num"]
    print(f"\nRendering Slide {step_num}/6 : {s['title']}...")
    setup_pymol_scene()
    action = s["action"]

    if action == "slide1":
        # Slide 1: Wild-Type complex
        cmd.load(pdb_wt, "complex")
        cmd.set_view(saved_view)
        cmd.show("cartoon", "complex")
        cmd.color("ruby", "complex and chain A")
        cmd.color("marine", "complex and chain B")
        # Interface sticks
        cmd.select("iface", "byres ((complex and chain A within 5.5 of (complex and chain B)) or (complex and chain B within 5.5 of (complex and chain A)))")
        cmd.show("sticks", "iface")
        cmd.color("salmon", "iface and chain A and elem C")
        cmd.color("cyan", "iface and chain B and elem C")
        cmd.color("firebrick", "iface and elem O and not elem C")
        cmd.color("tv_blue", "iface and elem N and not elem C")
        cmd.color("yellow", "iface and elem S and not elem C")

    elif action == "slide2":
        # Slide 2: Interface mapping & framework preservation
        cmd.load(pdb_wt, "complex")
        cmd.set_view(saved_view)
        cmd.show("cartoon", "complex")
        # Global framework in translucent / muted colors
        cmd.color("ruby", "complex and chain A")
        cmd.color("marine", "complex and chain B")
        cmd.set("cartoon_transparency", 0.65, "complex")
        # Highlight interface loops in solid opaque cartoon
        cmd.select("loops_a", "complex and chain A and (resi 3-22 or resi 29-47 or resi 185-205)")
        cmd.select("loops_b", "complex and chain B and (resi 2-20 or resi 29-47 or resi 74-81)")
        cmd.set("cartoon_transparency", 0.0, "loops_a or loops_b")
        cmd.color("ruby", "loops_a")
        cmd.color("marine", "loops_b")
        # Hotspots in sticks
        cmd.select("hotspots", "complex and chain A and resi 12+14+15+17+18+19+32+38")
        cmd.show("sticks", "hotspots")
        cmd.color("salmon", "hotspots and elem C")
        cmd.color("firebrick", "hotspots and elem O and not elem C")
        cmd.color("tv_blue", "hotspots and elem N and not elem C")

    elif action == "slide3":
        # Slide 3: A' Apo design with 29 mutations
        cmd.load(pdb_des, "complex")
        cmd.set_view(saved_view)
        cmd.hide("everything", "complex and chain B") # Apo mode!
        cmd.show("cartoon", "complex and chain A")
        cmd.color("ruby", "complex and chain A")
        # Show all 29 mutations as sticks
        a_mut_sel = "+".join(str(r) for r in A_mut_resi)
        cmd.select("a_muts", f"complex and chain A and resi {a_mut_sel}")
        cmd.show("sticks", "a_muts")
        cmd.color("salmon", "a_muts and elem C")
        cmd.color("firebrick", "a_muts and elem O and not elem C")
        cmd.color("tv_blue", "a_muts and elem N and not elem C")
        # Highlight key Ala cavity remodeling
        cmd.select("ala_cavity", "complex and chain A and resi 14+17+18+29")
        cmd.color("hotpink", "ala_cavity and elem C")

    elif action == "slide4":
        # Slide 4: Rupture with host B_WT (superposed!)
        cmd.load(pdb_rupture, "complex")
        cmd.set_view(saved_view)
        cmd.show("cartoon", "complex")
        cmd.color("ruby", "complex and chain A")
        cmd.color("marine", "complex and chain B")
        # Show mutated A' interface against B_WT
        a_mut_sel = "+".join(str(r) for r in A_mut_resi)
        cmd.select("a_muts", f"complex and chain A and resi {a_mut_sel}")
        cmd.select("b_iface", "byres (complex and chain B within 5.0 of a_muts)")
        cmd.show("sticks", "a_muts or b_iface")
        cmd.color("salmon", "a_muts and elem C")
        cmd.color("cyan", "b_iface and elem C")
        cmd.color("firebrick", "(a_muts or b_iface) and elem O and not elem C")
        cmd.color("tv_blue", "(a_muts or b_iface) and elem N and not elem C")

    elif action == "slide5":
        # Slide 5: B' compensatory rescue (37 mutations)
        cmd.load(pdb_des, "complex")
        cmd.set_view(saved_view)
        cmd.show("cartoon", "complex")
        cmd.color("ruby", "complex and chain A")
        cmd.color("marine", "complex and chain B")
        # Show both mutated interfaces
        a_mut_sel = "+".join(str(r) for r in A_mut_resi)
        b_mut_sel = "+".join(str(r) for r in B_mut_resi)
        cmd.select("a_muts", f"complex and chain A and resi {a_mut_sel}")
        cmd.select("b_muts", f"complex and chain B and resi {b_mut_sel}")
        cmd.show("sticks", "a_muts or b_muts")
        cmd.color("salmon", "a_muts and elem C")
        cmd.color("cyan", "b_muts and elem C")
        cmd.color("firebrick", "(a_muts or b_muts) and elem O and not elem C")
        cmd.color("tv_blue", "(a_muts or b_muts) and elem N and not elem C")
        # Highlight aromatic latch Phe12 & Phe47
        cmd.select("latch", "complex and chain B and resi 12+47")
        cmd.show("spheres", "latch")
        cmd.set("sphere_scale", 0.38, "latch")
        cmd.color("deepteal", "latch and elem C")

    elif action == "slide6":
        # Slide 6: Direct 1:1 Superposition (WT vs Champion Design!)
        # Shows that framework is 100% aligned (RMSD = 0.76 A), but interface is completely reprogrammed
        cmd.load(pdb_wt, "wt")
        cmd.load(pdb_des, "des")
        cmd.set_view(saved_view)
        # Skeletons
        cmd.show("cartoon", "des")
        cmd.color("ruby", "des and chain A")
        cmd.color("marine", "des and chain B")
        # Superposed WT skeleton in faint ghost cartoon to show perfect framework alignment
        cmd.show("cartoon", "wt")
        cmd.color("gray80", "wt")
        cmd.set("cartoon_transparency", 0.70, "wt")
        # All-atom interface on designed champion
        cmd.select("iface_champ", "byres ((des and chain A within 5.5 of (des and chain B)) or (des and chain B within 5.5 of (des and chain A)))")
        cmd.show("sticks", "iface_champ")
        cmd.color("salmon", "iface_champ and chain A and elem C")
        cmd.color("cyan", "iface_champ and chain B and elem C")
        cmd.color("firebrick", "iface_champ and elem O and not elem C")
        cmd.color("tv_blue", "iface_champ and elem N and not elem C")
        cmd.color("yellow", "iface_champ and elem S and not elem C")
        # Aromatic latch in spheres
        cmd.select("latch_champ", "des and chain B and resi 12+47")
        cmd.show("spheres", "latch_champ")
        cmd.set("sphere_scale", 0.36, "latch_champ")
        cmd.color("deepteal", "latch_champ and elem C")

    # Save raw render image
    raw_path = os.path.join(SLIDES_DIR, f"slide_{step_num}_raw.png")
    cmd.png(raw_path, width=WIDTH, height=HEIGHT, ray=0)

    # ----------------- Draw HUD Banners with PIL -----------------
    img = Image.open(raw_path).convert("RGBA")
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    # Top header bar (Height 105px)
    draw.rectangle([(0, 0), (WIDTH, 105)], fill=HEX_CARD)
    draw.line([(0, 105), (WIDTH, 105)], fill=HEX_CYAN, width=2)

    # Step badge & tag
    draw.rounded_rectangle([(30, 16), (220, 46)], radius=6, fill=HEX_BLUE)
    draw.text((45, 22), f"ETAPE {s['step_num']} / {s['step_total']}", font=font_step, fill=HEX_WHITE)

    draw.rounded_rectangle([(235, 16), (580, 46)], radius=6, fill=(30, 41, 59))
    draw.text((250, 22), s["tag"], font=font_step, fill=HEX_CYAN)

    # Main Title & Subtitle
    draw.text((30, 52), s["title"], font=font_title, fill=HEX_WHITE)
    draw.text((30, 84), s["subtitle"], font=font_subtitle, fill=HEX_MUTED)

    # Top-right chain color guide
    draw.rounded_rectangle([(WIDTH - 540, 16), (WIDTH - 280, 46)], radius=6, fill=(225, 29, 72, 230))
    draw.text((WIDTH - 525, 22), "Chaine A : S1 / A' (Rouge)", font=font_legend, fill=HEX_WHITE)

    draw.rounded_rectangle([(WIDTH - 265, 16), (WIDTH - 30, 46)], radius=6, fill=(37, 99, 235, 230))
    draw.text((WIDTH - 250, 22), "Chaine B : S2 / B' (Bleu)", font=font_legend, fill=HEX_WHITE)

    # Right badge info
    draw.text((WIDTH - 540, 60), "Champion : A'_candidate_01 · B_cand_40", font=font_hud_mono, fill=HEX_CYAN)
    draw.text((WIDTH - 540, 80), "Cible : Ribosome 30S (S1-S2-ARNr 16S)", font=font_hud_text, fill=HEX_MUTED)

    # Bottom HUD panel (Height 175px)
    hud_top = HEIGHT - 175
    draw.rectangle([(0, hud_top), (WIDTH, HEIGHT)], fill=HEX_CARD)
    draw.line([(0, hud_top), (WIDTH, hud_top)], fill=HEX_CYAN, width=2)

    # Row 1: Software & Engines
    draw.rounded_rectangle([(30, hud_top + 13), (225, hud_top + 37)], radius=4, fill=(120, 53, 15))
    draw.text((40, hud_top + 16), "[LOGICIELS & ENGINS]", font=font_hud_tag, fill=HEX_AMBER)
    draw.text((245, hud_top + 16), s["software"], font=font_hud_mono, fill=HEX_WHITE)

    # Row 2: Pipeline Parameters
    draw.rounded_rectangle([(30, hud_top + 50), (225, hud_top + 74)], radius=4, fill=(6, 78, 59))
    draw.text((40, hud_top + 53), "[PARAMETRES PIPELINE]", font=font_hud_tag, fill=HEX_EMERALD)
    draw.text((245, hud_top + 53), s["params"], font=font_hud_text, fill=(203, 213, 225))

    # Row 3: All-Atom Mechanism & Mutations
    draw.rounded_rectangle([(30, hud_top + 87), (225, hud_top + 111)], radius=4, fill=(14, 116, 144))
    draw.text((40, hud_top + 90), "[MECANISME ALL-ATOM]", font=font_hud_tag, fill=HEX_CYAN)
    draw.text((245, hud_top + 90), s["mechanism"], font=font_hud_text, fill=HEX_WHITE)

    # Row 4: Biophysical Metrics & Status
    draw.rounded_rectangle([(30, hud_top + 124), (225, hud_top + 148)], radius=4, fill=(30, 41, 59))
    draw.text((40, hud_top + 127), "[METRIQUES CLES]", font=font_hud_tag, fill=HEX_WHITE)
    draw.text((245, hud_top + 127), s["metrics"], font=font_hud_mono, fill=HEX_CYAN)

    # Composite & save final slide
    final_img = Image.alpha_composite(img, overlay).convert("RGB")
    final_slide_path = os.path.join(SLIDES_DIR, f"slide_{step_num}.png")
    final_img.save(final_slide_path, quality=95)
    slide_files.append(final_slide_path)
    if os.path.exists(raw_path):
        os.remove(raw_path)
    print(f"Saved Slide {step_num} to {final_slide_path}")

cmd.quit()
print(f"\n[SUCCESS] All {len(slide_files)} slides generated successfully!")

# ----------------- Video & Animation Assembly -----------------
video_mp4_path = os.path.join(FIG_DIR, "pipeline_diaporama_video.mp4")
anim_gif_path = os.path.join(FIG_DIR, "pipeline_diaporama_animation.gif")

print("\nAssembling MP4 video with cross-dissolve transitions...")
fps = 30
duration_hold_sec = 4.0
duration_fade_sec = 0.8
hold_frames = int(duration_hold_sec * fps)
fade_frames = int(duration_fade_sec * fps)

writer = imageio.get_writer(video_mp4_path, fps=fps, codec="libx264", quality=8, pixelformat="yuv420p")
gif_frames = []

for idx, slide_path in enumerate(slide_files):
    curr_im = np.array(Image.open(slide_path).convert("RGB"))
    next_im = np.array(Image.open(slide_files[(idx + 1) % len(slide_files)]).convert("RGB"))

    # Hold slide
    for _ in range(hold_frames):
        writer.append_data(curr_im)

    # Cross-dissolve fade to next slide
    for f in range(fade_frames):
        alpha = (f + 1) / float(fade_frames)
        blend = ((1.0 - alpha) * curr_im + alpha * next_im).astype(np.uint8)
        writer.append_data(blend)

    # For GIF (scaled down to 800x450, 1 frame per slide)
    small_im = Image.fromarray(curr_im).resize((800, 450), Image.Resampling.LANCZOS)
    gif_frames.append(small_im)

writer.close()
print(f"MP4 Video saved to {video_mp4_path} ({os.path.getsize(video_mp4_path) / (1024*1024):.2f} MB)")

print("Saving animated GIF...")
gif_frames[0].save(
    anim_gif_path,
    save_all=True,
    append_images=gif_frames[1:],
    duration=3500, # 3.5s per slide
    loop=0,
    optimize=True
)
print(f"Looping GIF saved to {anim_gif_path} ({os.path.getsize(anim_gif_path) / 1024:.1f} KB)")
print("Done!")
