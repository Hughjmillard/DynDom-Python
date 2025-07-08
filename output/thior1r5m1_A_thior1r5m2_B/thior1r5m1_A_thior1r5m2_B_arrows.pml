# DynDom Arrow Visualization Script
# Load main protein structure
load output\thior1r5m1_A_thior1r5m2_B\thior1r5m1_A_thior1r5m2_B.pdb

# Load arrows
load thior1r5m1_A_thior1r5m2_B_arrows.pdb

# Basic protein display
hide everything, output\thior1r5m1_A_thior1r5m2_B\thior1r5m1_A_thior1r5m2_B
show cartoon, output\thior1r5m1_A_thior1r5m2_B\thior1r5m1_A_thior1r5m2_B
color gray80, output\thior1r5m1_A_thior1r5m2_B\thior1r5m1_A_thior1r5m2_B

# Hide arrow atoms initially
hide everything, thior1r5m1_A_thior1r5m2_B_arrows

# Arrow 1: Domain 0 (moving) relative to Domain 1 (fixed)
# Shaft color: blue (fixed domain), Head color: red (moving domain)
# Rotation: 66.4°, Translation: -1.5Å

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
bg_color white
set depth_cue, 0
set ray_shadows, 1
set ray_shadow_decay_factor, 0.1

# Better lighting for 3D arrow heads
set ambient, 0.2
set direct, 0.8
set reflect, 0.5
set shininess, 10

# Center view
zoom all
orient

# Clean up selections
delete shaft_*
delete head_*

print 'DynDom arrows with 3D heads loaded successfully!'
print 'Fixed domain: 1 (red)'
print 'Moving domain 0: Chain A, red shaft (fixed) with blue head (moving), 66.4° rotation'