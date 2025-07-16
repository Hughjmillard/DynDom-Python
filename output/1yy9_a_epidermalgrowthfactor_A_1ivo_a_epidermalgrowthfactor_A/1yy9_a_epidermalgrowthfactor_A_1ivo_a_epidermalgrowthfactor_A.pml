# DynDom Complete Visualization
# Domain-colored structure with screw axis arrows
reinitialize
load 1yy9_a_epidermalgrowthfactor_A_1ivo_a_epidermalgrowthfactor_A.pdb
bg_color white
color grey

# === DOMAIN STRUCTURE COLORING ===
select fixed_domain, resi 6-6
select fixed_domain, fixed_domain + resi 229-229
select fixed_domain, fixed_domain + resi 231-231
select fixed_domain, fixed_domain + resi 234-272
select fixed_domain, fixed_domain + resi 276-279
color blue, fixed_domain

select moving_domain_0, resi 311-507
color red, moving_domain_0

select moving_domain_1, resi 3-4
select moving_domain_1, moving_domain_1 + resi 8-227
color yellow, moving_domain_1

select moving_domain_2, resi 274-274
select moving_domain_2, moving_domain_2 + resi 281-308
color pink, moving_domain_2

# Color bending residues
select bending_residues_1, resi 310-310
color green, bending_residues_1
select bending_residues_2, resi 508-508
color green, bending_residues_2
select bending_residues_3, resi 2-2
color green, bending_residues_3
select bending_residues_4, resi 5-5
color green, bending_residues_4
select bending_residues_5, resi 7-7
color green, bending_residues_5
select bending_residues_6, resi 228-228
color green, bending_residues_6
select bending_residues_7, resi 230-230
color green, bending_residues_7
select bending_residues_8, resi 232-233
color green, bending_residues_8
select bending_residues_9, resi 273-273
color green, bending_residues_9
select bending_residues_10, resi 275-275
color green, bending_residues_10
select bending_residues_11, resi 280-280
color green, bending_residues_11
select bending_residues_12, resi 309-309
color green, bending_residues_12

set dash_gap, 0
set dash_radius, 0.2

# === SCREW AXIS ARROWS ===
load output
load 1yy9_a_epidermalgrowthfactor_A_1ivo_a_epidermalgrowthfactor_A_arrows.pdb

# Basic protein display
hide everything, output
show cartoon, output
color gray80, output

# Hide arrow atoms initially
hide everything, 1yy9_a_epidermalgrowthfactor_A_1ivo_a_epidermalgrowthfactor_A_arrows

# Arrow 1: Domain 0 (moving) relative to Domain 3 (fixed)
# Shaft color: blue (fixed domain), Head color: red (moving domain)
# Rotation: 136.2°

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

# Arrow 2: Domain 1 (moving) relative to Domain 3 (fixed)
# Shaft color: blue (fixed domain), Head color: yellow (moving domain)
# Rotation: 138.9°

# Select shaft and head atoms by chain and residue
select shaft_2, chain B and resn SHF and resi 150
select head_2, chain B and resn ARH and resi 170

# Display shaft as thick licorice stick (FIXED domain color: blue)
show sticks, shaft_2
color blue, shaft_2
set stick_radius, 0.3, shaft_2

# Display arrow head as clean cone (MOVING domain color: yellow)
show sticks, head_2
color yellow, head_2
set stick_radius, 0.25, head_2

# Connect atoms ONLY within each section
bond shaft_2, shaft_2
bond head_2, head_2

# Arrow 3: Domain 2 (moving) relative to Domain 3 (fixed)
# Shaft color: blue (fixed domain), Head color: pink (moving domain)
# Rotation: 124.3°

# Select shaft and head atoms by chain and residue
select shaft_3, chain C and resn SHF and resi 200
select head_3, chain C and resn ARH and resi 220

# Display shaft as thick licorice stick (FIXED domain color: blue)
show sticks, shaft_3
color blue, shaft_3
set stick_radius, 0.3, shaft_3

# Display arrow head as clean cone (MOVING domain color: pink)
show sticks, head_3
color pink, head_3
set stick_radius, 0.25, head_3

# Connect atoms ONLY within each section
bond shaft_3, shaft_3
bond head_3, head_3

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
print 'Fixed domain: 3 (blue)'
print 'Moving domain 0: red, rotation 136.2°'
print 'Moving domain 1: yellow, rotation 138.9°'
print 'Moving domain 2: pink, rotation 124.3°'
