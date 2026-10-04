# PyMOL Visualization Script - OrthoPPInterface Top Champions
reinitialize
bg_color white
set ray_shadows, 0
set cartoon_transparency, 0.1
set stick_radius, 0.25

load ../structures_pdb/01_A_prime_candidate_01_B_cand_40_design.pdb, rank1_A_prime_candidate_01_B_ca
color cyan, rank1_A_prime_candidate_01_B_ca and chain A
color sandybrown, rank1_A_prime_candidate_01_B_ca and chain B
select rank1_A_prime_candidate_01_B_ca_mutA, rank1_A_prime_candidate_01_B_ca and chain A and resi 3+4+12+13+14+15+16+17+18+19+29+32+33+34+37+38+39+40+41+70+71+72+185+186+201+202+203+204+205
show sticks, rank1_A_prime_candidate_01_B_ca_mutA
color deepteal, rank1_A_prime_candidate_01_B_ca_mutA
select rank1_A_prime_candidate_01_B_ca_mutB, rank1_A_prime_candidate_01_B_ca and chain B and resi 2+3+4+5+6+7+8+9+10+11+12+16+17+18+19+20+29+30+32+35+37+39+40+41+42+43+44+45+47+74+75+76+77+78+79+80+81
show sticks, rank1_A_prime_candidate_01_B_ca_mutB
color firebrick, rank1_A_prime_candidate_01_B_ca_mutB
disable rank1_A_prime_candidate_01_B_ca
load ../structures_pdb/02_A_prime_candidate_17_B_cand_29_design.pdb, rank2_A_prime_candidate_17_B_ca
color cyan, rank2_A_prime_candidate_17_B_ca and chain A
color sandybrown, rank2_A_prime_candidate_17_B_ca and chain B
select rank2_A_prime_candidate_17_B_ca_mutA, rank2_A_prime_candidate_17_B_ca and chain A and resi 3+4+10+13+15+16+17+18+19+29+32+34+37+38+39+41+70+71+75+185+186+201+202+203+204+205
show sticks, rank2_A_prime_candidate_17_B_ca_mutA
color deepteal, rank2_A_prime_candidate_17_B_ca_mutA
select rank2_A_prime_candidate_17_B_ca_mutB, rank2_A_prime_candidate_17_B_ca and chain B and resi 2+3+4+5+6+7+8+10+11+12+16+17+18+19+29+30+32+35+37+39+40+41+42+43+44+45+47+74+75+76+77+78+79
show sticks, rank2_A_prime_candidate_17_B_ca_mutB
color firebrick, rank2_A_prime_candidate_17_B_ca_mutB
disable rank2_A_prime_candidate_17_B_ca
load ../structures_pdb/03_A_prime_candidate_01_B_cand_36_design.pdb, rank3_A_prime_candidate_01_B_ca
color cyan, rank3_A_prime_candidate_01_B_ca and chain A
color sandybrown, rank3_A_prime_candidate_01_B_ca and chain B
select rank3_A_prime_candidate_01_B_ca_mutA, rank3_A_prime_candidate_01_B_ca and chain A and resi 3+4+12+13+14+15+16+17+18+19+29+32+33+34+37+38+39+40+41+70+71+72+185+186+201+202+203+204+205
show sticks, rank3_A_prime_candidate_01_B_ca_mutA
color deepteal, rank3_A_prime_candidate_01_B_ca_mutA
select rank3_A_prime_candidate_01_B_ca_mutB, rank3_A_prime_candidate_01_B_ca and chain B and resi 2+3+4+5+6+7+8+9+10+12+16+17+18+19+20+29+30+32+35+37+39+40+42+43+44+45+47+74+75+76+78+79+80+81
show sticks, rank3_A_prime_candidate_01_B_ca_mutB
color firebrick, rank3_A_prime_candidate_01_B_ca_mutB
disable rank3_A_prime_candidate_01_B_ca
load ../structures_pdb/04_A_prime_candidate_13_B_cand_22_design.pdb, rank4_A_prime_candidate_13_B_ca
color cyan, rank4_A_prime_candidate_13_B_ca and chain A
color sandybrown, rank4_A_prime_candidate_13_B_ca and chain B
select rank4_A_prime_candidate_13_B_ca_mutA, rank4_A_prime_candidate_13_B_ca and chain A and resi 3+4+12+14+15+16+17+18+19+29+32+33+34+35+38+39+41+70+71+72+75+185+186+201+202+203+204+205
show sticks, rank4_A_prime_candidate_13_B_ca_mutA
color deepteal, rank4_A_prime_candidate_13_B_ca_mutA
select rank4_A_prime_candidate_13_B_ca_mutB, rank4_A_prime_candidate_13_B_ca and chain B and resi 1+2+3+4+5+6+7+8+9+10+12+13+16+17+18+19+20+30+32+39+40+42+43+44+45+75+76+77+78+79+80+81
show sticks, rank4_A_prime_candidate_13_B_ca_mutB
color firebrick, rank4_A_prime_candidate_13_B_ca_mutB
disable rank4_A_prime_candidate_13_B_ca
load ../structures_pdb/05_A_prime_candidate_11_B_cand_33_design.pdb, rank5_A_prime_candidate_11_B_ca
color cyan, rank5_A_prime_candidate_11_B_ca and chain A
color sandybrown, rank5_A_prime_candidate_11_B_ca and chain B
select rank5_A_prime_candidate_11_B_ca_mutA, rank5_A_prime_candidate_11_B_ca and chain A and resi 3+4+10+13+15+16+17+18+19+29+32+34+35+37+38+39+40+41+70+71+75+185+186+201+202+203+204+205
show sticks, rank5_A_prime_candidate_11_B_ca_mutA
color deepteal, rank5_A_prime_candidate_11_B_ca_mutA
select rank5_A_prime_candidate_11_B_ca_mutB, rank5_A_prime_candidate_11_B_ca and chain B and resi 1+2+4+5+6+7+9+10+12+16+17+18+19+20+29+32+35+37+39+40+41+42+43+44+45+47+74+75+76+77+78+79+80
show sticks, rank5_A_prime_candidate_11_B_ca_mutB
color firebrick, rank5_A_prime_candidate_11_B_ca_mutB
disable rank5_A_prime_candidate_11_B_ca

enable rank1_A_prime_candidate_01_B_ca
zoom rank1_A_prime_candidate_01_B_ca
orient
