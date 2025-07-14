# DynDom Complete Visualization
# Domain-colored structure with screw axis arrows
reinitialize
load 4ake_A_2eck_B.pdb
bg_color white
color grey

# === DOMAIN STRUCTURE COLORING ===
select fixed_domain, resi 1-25
select fixed_domain, fixed_domain + resi 63-111
select fixed_domain, fixed_domain + resi 169-210
color blue, fixed_domain

select moving_domain_0, resi 116-152
color red, moving_domain_0

select moving_domain_1, resi 29-58
color yellow, moving_domain_1

# Color bending residues
select bending_residues_1, resi 112-115
color green, bending_residues_1
select bending_residues_2, resi 153-168
color green, bending_residues_2
select bending_residues_3, resi 26-28
color green, bending_residues_3
select bending_residues_4, resi 59-62
color green, bending_residues_4

set dash_gap, 0
set dash_radius, 0.2

# === SCREW AXIS ARROWS ===
load 4ake_A_2eck_B_arrows.pdb

# Hide arrow atoms initially
hide everything, 4ake_A_2eck_B_arrows

# Arrow 1: Domain 0 (moving) relative to Domain 2 (fixed)
# Shaft color: blue (fixed domain), Head color: red (moving domain)
# Rotation: 53.0°

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

# Arrow 2: Domain 1 (moving) relative to Domain 2 (fixed)
# Shaft color: blue (fixed domain), Head color: yellow (moving domain)
# Rotation: 46.1°

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

# Disable automatic bonding between different chains
set auto_bond, 0

# Make arrows more prominent
set stick_transparency, 0.0
set stick_quality, 15
set sphere_quality, 3
set surface_quality, 2

# Better lighting for 3D arrow heads
set ambient, 0.2
set direct, 0.8
set reflect, 0.5
set shininess, 10

# Clean up selections
delete shaft_*
delete head_*

# === FINAL SETTINGS ===
show cartoon
set cartoon_transparency, 0.1
set stick_transparency, 0.0
set stick_quality, 15
zoom all
orient

# Cleanup selections
delete fixed_domain
delete moving_domain_*
delete bending_residues
delete arrow_*

print 'DynDom complete visualization loaded!'
print 'Fixed domain: 2 (blue)'
print 'Moving domain 0: red, rotation 53.0°'
print 'Moving domain 1: yellow, rotation 46.1°'
