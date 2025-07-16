# DynDom Complete Visualization
# Domain-colored structure with screw axis arrows
reinitialize
load 2dcy_e_2qz3_a_A_2dcy_e_2qz3_a_B.pdb
bg_color white
color grey

# === DOMAIN STRUCTURE COLORING ===
select fixed_domain, resi 2-40
select fixed_domain, fixed_domain + resi 48-64
select fixed_domain, fixed_domain + resi 80-91
select fixed_domain, fixed_domain + resi 129-141
select fixed_domain, fixed_domain + resi 150-156
select fixed_domain, fixed_domain + resi 166-181
color blue, fixed_domain

select moving_domain_0, resi 68-77
select moving_domain_0, moving_domain_0 + resi 94-97
select moving_domain_0, moving_domain_0 + resi 102-126
select moving_domain_0, moving_domain_0 + resi 144-147
select moving_domain_0, moving_domain_0 + resi 159-163
color red, moving_domain_0

# Color bending residues
select bending_residues_1, resi 41-47
color green, bending_residues_1
select bending_residues_2, resi 65-67
color green, bending_residues_2
select bending_residues_3, resi 78-79
color green, bending_residues_3
select bending_residues_4, resi 92-93
color green, bending_residues_4
select bending_residues_5, resi 98-101
color green, bending_residues_5
select bending_residues_6, resi 127-128
color green, bending_residues_6
select bending_residues_7, resi 142-143
color green, bending_residues_7
select bending_residues_8, resi 148-149
color green, bending_residues_8
select bending_residues_9, resi 157-158
color green, bending_residues_9
select bending_residues_10, resi 164-165
color green, bending_residues_10

set dash_gap, 0
set dash_radius, 0.2

# === SCREW AXIS ARROWS ===
load output
load 2dcy_e_2qz3_a_A_2dcy_e_2qz3_a_B_arrows.pdb

# Basic protein display
hide everything, output
show cartoon, output
color gray80, output

# Hide arrow atoms initially
hide everything, 2dcy_e_2qz3_a_A_2dcy_e_2qz3_a_B_arrows

# Arrow 1: Domain 0 (moving) relative to Domain 1 (fixed)
# Shaft color: blue (fixed domain), Head color: red (moving domain)
# Rotation: 6.8°

# Select shaft and head atoms by chain and residue
select shaft_1, chain A and resn SHF and resi 100
select head_1, chain A and resn ARH and resi 120

# Display shaft as thick licorice stick (FIXED domain color: blue)
show sticks, shaft_1
color blue, shaft_1
set stick_radius, 0.3, shaft_1

# Display arrow head as clean cone (MOVING domain color: red)
show sticks, head_1
color red, head_1
set stick_radius, 0.25, head_1

# Connect atoms ONLY within each section
bond shaft_1, shaft_1
bond head_1, head_1

# Disable automatic bonding between different chains
set auto_bond, 0

# Make arrows more prominent
set stick_transparency, 0.0
set stick_quality, 15
set sphere_quality, 3
set surface_quality, 2

# Final settings
set depth_cue, 0
set ray_shadows, 1
set ray_shadow_decay_factor, 0.1

# Better lighting for 3D arrow heads
set ambient, 0.2
set direct, 0.8
set reflect, 0.5
set shininess, 10

# Clean up selections
delete shaft_*
delete head_*

# === FINAL SETTINGS ===
set stick_transparency, 0.0
set stick_quality, 15
zoom all
orient

# Cleanup selections
delete bending_residues
delete arrow_*

print 'DynDom complete visualization loaded!'
print 'Fixed domain: 1 (blue)'
print 'Moving domain 0: red, rotation 6.8°'
