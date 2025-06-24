# DynDom Arrow Visualization Script
# Load main protein structure
load output\1adg_A_2ohx_B\1adg_A_2ohx_B.pdb

# Load arrows
load 1adg_A_2ohx_B_arrows.pdb

# Basic protein display
hide everything, output\1adg_A_2ohx_B\1adg_A_2ohx_B
show cartoon, output\1adg_A_2ohx_B\1adg_A_2ohx_B
color gray80, output\1adg_A_2ohx_B\1adg_A_2ohx_B

# Hide arrow atoms initially
hide everything, 1adg_A_2ohx_B_arrows

# Arrow 1: Domain 1 -> Domain 0
# Rotation: 10.0°, Translation: -0.3Å

# Select shaft and head atoms
select shaft_1, resn SHF and resi 10
select head_1, resn ARH and resi 11

# Display shaft (fixed domain color: blue)
show spheres, shaft_1
show sticks, shaft_1
color blue, shaft_1
set sphere_scale, 0.8, shaft_1
set stick_radius, 0.3, shaft_1

# Display head (moving domain color: red)
show spheres, head_1
show sticks, head_1
color red, head_1
set sphere_scale, 1.2, head_1
set stick_radius, 0.5, head_1

# Connect atoms within each section
bond shaft_1, shaft_1
bond head_1, head_1
# Connect shaft to head
bond (resi 10 and name CA), (resi 11 and name CA)

# Make arrows more prominent
set stick_transparency, 0.0
set sphere_transparency, 0.0
set stick_quality, 15
set sphere_quality, 3

# Final settings
bg_color white
set depth_cue, 0
set ray_shadows, 0

# Show domain assignments on protein
color blue, output\1adg_A_2ohx_B\1adg_A_2ohx_B and not (resn SHF or resn ARH)

# Center view
zoom all
orient

# Clean up selections
delete shaft_*
delete head_*

print 'DynDom arrows loaded successfully!'
print 'Fixed domain: 0 (blue)'
print 'Moving domain 1: red arrow, 10.0° rotation'