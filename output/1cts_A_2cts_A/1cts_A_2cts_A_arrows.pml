# DynDom Arrow Visualization Script
# Load main protein structure
load output\1cts_A_2cts_A\1cts_A_2cts_A.pdb

# Load arrows
load 1cts_A_2cts_A_arrows.pdb

# Basic protein display
hide everything, output\1cts_A_2cts_A\1cts_A_2cts_A
show cartoon, output\1cts_A_2cts_A\1cts_A_2cts_A
color gray80, output\1cts_A_2cts_A\1cts_A_2cts_A

# Hide arrow atoms initially
hide everything, 1cts_A_2cts_A_arrows

# Arrow 1: Domain 0 -> Domain 1
select shaft_1, resn SHF and resi 1000-1050
select head_1, resn ARH and resi 1051-1099

# Display shaft (fixed domain color)
show spheres, shaft_1
show sticks, shaft_1
color red, shaft_1
set sphere_scale, 0.3, shaft_1
set stick_radius, 0.1, shaft_1

# Display head (moving domain color)
show spheres, head_1
show sticks, head_1
color blue, head_1
set sphere_scale, 0.5, head_1
set stick_radius, 0.2, head_1

# Connect atoms
bond shaft_1, shaft_1
bond head_1, head_1

# Final settings
bg_color white
zoom all

# Clean up selections
delete shaft_*
delete head_*

print 'DynDom arrows loaded successfully!'