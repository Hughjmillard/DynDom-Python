# DynDom Complete Visualization
# Domain-colored structure with screw axis arrows
reinitialize
load 1adg_A_2ohx_B.pdb
bg_color white
color grey

# === DOMAIN STRUCTURE COLORING ===
select fixed_domain, resi 1-170
select fixed_domain, fixed_domain + resi 317-368
color blue, fixed_domain

select moving_domain_1, resi 175-282
select moving_domain_1, moving_domain_1 + resi 295-311
color red, moving_domain_1

# Color bending residues
select bending_residues_1, resi 171-174
color green, bending_residues_1
select bending_residues_2, resi 283-294
color green, bending_residues_2
select bending_residues_3, resi 312-316
color green, bending_residues_3

set dash_gap, 0
set dash_radius, 0.2

# === SCREW AXIS ARROWS ===
load output
load 1adg_A_2ohx_B_arrows.pdb

# Basic protein display
hide everything, output
show cartoon, output
color gray80, output

# Hide arrow atoms initially
hide everything, 1adg_A_2ohx_B_arrows

# Arrow 1: Domain 1 (moving) relative to Domain 0 (fixed)
# Shaft color: blue (fixed domain), Head color: red (moving domain)
# Rotation: 10.0°

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
print 'Fixed domain: 0 (blue)'
print 'Moving domain 1: red, rotation 10.0°'
