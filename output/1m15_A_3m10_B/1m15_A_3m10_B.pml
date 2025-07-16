# DynDom Complete Visualization
# Domain-colored structure with screw axis arrows
reinitialize
load 1m15_A_3m10_B.pdb
bg_color white
color grey

# === DOMAIN STRUCTURE COLORING ===
select fixed_domain, resi 2-94
select fixed_domain, fixed_domain + resi 125-151
select fixed_domain, fixed_domain + resi 163-167
select fixed_domain, fixed_domain + resi 190-205
select fixed_domain, fixed_domain + resi 218-226
select fixed_domain, fixed_domain + resi 257-273
color blue, fixed_domain

select moving_domain_1, resi 98-122
select moving_domain_1, moving_domain_1 + resi 156-160
select moving_domain_1, moving_domain_1 + resi 211-211
select moving_domain_1, moving_domain_1 + resi 231-251
select moving_domain_1, moving_domain_1 + resi 277-310
select moving_domain_1, moving_domain_1 + resi 321-353
color red, moving_domain_1

select moving_domain_2, resi 172-181
select moving_domain_2, moving_domain_2 + resi 212-215
color yellow, moving_domain_2

# Color bending residues
select bending_residues_1, resi 95-97
color green, bending_residues_1
select bending_residues_2, resi 123-124
color green, bending_residues_2
select bending_residues_3, resi 152-155
color green, bending_residues_3
select bending_residues_4, resi 161-162
color green, bending_residues_4
select bending_residues_5, resi 227-230
color green, bending_residues_5
select bending_residues_6, resi 252-256
color green, bending_residues_6
select bending_residues_7, resi 274-276
color green, bending_residues_7
select bending_residues_8, resi 168-171
color green, bending_residues_8
select bending_residues_9, resi 182-189
color green, bending_residues_9
select bending_residues_10, resi 206-210
color green, bending_residues_10
select bending_residues_11, resi 216-217
color green, bending_residues_11

set dash_gap, 0
set dash_radius, 0.2

# === SCREW AXIS ARROWS ===
load output
load 1m15_A_3m10_B_arrows.pdb

# Basic protein display
hide everything, output
show cartoon, output
color gray80, output

# Hide arrow atoms initially
hide everything, 1m15_A_3m10_B_arrows

# Arrow 1: Domain 1 (moving) relative to Domain 0 (fixed)
# Shaft color: blue (fixed domain), Head color: red (moving domain)
# Rotation: 20.4°

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

# Arrow 2: Domain 2 (moving) relative to Domain 0 (fixed)
# Shaft color: blue (fixed domain), Head color: yellow (moving domain)
# Rotation: 20.9°

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
print 'Moving domain 1: red, rotation 20.4°'
print 'Moving domain 2: yellow, rotation 20.9°'
