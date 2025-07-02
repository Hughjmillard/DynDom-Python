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

# Arrow 1: Domain 0 (moving) relative to Domain 2 (fixed)
# Shaft color: yellow (fixed domain), Head color: blue (moving domain)
# Rotation: 53.0°, Translation: 1.2Å

# Select shaft and head atoms by chain and residue
select shaft_1, chain A and resn SHF and resi 100
select head_1, chain A and resn ARH and resi 120

# Display shaft as thick licorice stick (FIXED domain color: yellow)
show sticks, shaft_1
color yellow, shaft_1
set stick_radius, 0.3, shaft_1

# Display arrow head as clean cone (MOVING domain color: blue)
show sticks, head_1
color blue, head_1
set stick_radius, 0.25, head_1

# Connect atoms ONLY within each section
bond shaft_1, shaft_1
bond head_1, head_1

# Arrow 2: Domain 1 (moving) relative to Domain 2 (fixed)
# Shaft color: yellow (fixed domain), Head color: red (moving domain)
# Rotation: 46.1°, Translation: 1.4Å

# Select shaft and head atoms by chain and residue
select shaft_2, chain B and resn SHF and resi 150
select head_2, chain B and resn ARH and resi 170

# Display shaft as thick licorice stick (FIXED domain color: yellow)
show sticks, shaft_2
color yellow, shaft_2
set stick_radius, 0.3, shaft_2

# Display arrow head as clean cone (MOVING domain color: red)
show sticks, head_2
color red, head_2
set stick_radius, 0.25, head_2

# Connect atoms ONLY within each section
bond shaft_2, shaft_2
bond head_2, head_2

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
print 'Fixed domain: 2 (yellow)'
print 'Moving domain 0: Chain A, yellow shaft (fixed) with blue head (moving), 53.0° rotation'
print 'Moving domain 1: Chain B, yellow shaft (fixed) with red head (moving), 46.1° rotation'