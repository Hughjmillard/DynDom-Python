# DynDom Arrow Visualization Script
# Load main protein structure
load output\4ake_A_2eck_B\4ake_A_2eck_B.pdb

# Load arrows
load 4ake_A_2eck_B_arrows.pdb

# Basic protein display
hide everything, output\4ake_A_2eck_B\4ake_A_2eck_B
show cartoon, output\4ake_A_2eck_B\4ake_A_2eck_B
color gray80, output\4ake_A_2eck_B\4ake_A_2eck_B

# Hide arrow atoms initially
hide everything, 4ake_A_2eck_B_arrows

# Arrow 1: Domain 0 -> Domain 2
select shaft_1, resn SHF and resi 1000-1050
select head_1, resn ARH and resi 1051-1099

# Display shaft (fixed domain color)
show spheres, shaft_1
show sticks, shaft_1
color yellow, shaft_1
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

# Arrow 2: Domain 1 -> Domain 2
select shaft_2, resn SHF and resi 1100-1150
select head_2, resn ARH and resi 1151-1199

# Display shaft (fixed domain color)
show spheres, shaft_2
show sticks, shaft_2
color yellow, shaft_2
set sphere_scale, 0.3, shaft_2
set stick_radius, 0.1, shaft_2

# Display head (moving domain color)
show spheres, head_2
show sticks, head_2
color red, head_2
set sphere_scale, 0.5, head_2
set stick_radius, 0.2, head_2

# Connect atoms
bond shaft_2, shaft_2
bond head_2, head_2

# Final settings
bg_color white
zoom all

# Clean up selections
delete shaft_*
delete head_*

print 'DynDom arrows loaded successfully!'