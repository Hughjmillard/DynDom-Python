# DynDom Complete Visualization
# Domain-colored structure with screw axis arrows
reinitialize
load 1yy9_a_epidermalgrowthfactor_A_1ivo_a_epidermalgrowthfactor_A.pdb
bg_color white
color grey

# === DOMAIN STRUCTURE COLORING ===
select fixed_domain, resi 274-274
select fixed_domain, fixed_domain + resi 281-305
color blue, fixed_domain

select moving_domain_0, resi 311-508
color red, moving_domain_0

select moving_domain_1, resi 2-5
select moving_domain_1, moving_domain_1 + resi 7-228
select moving_domain_1, moving_domain_1 + resi 230-230
select moving_domain_1, moving_domain_1 + resi 232-233
color yellow, moving_domain_1

select moving_domain_3, resi 6-6
select moving_domain_3, moving_domain_3 + resi 229-229
select moving_domain_3, moving_domain_3 + resi 231-231
select moving_domain_3, moving_domain_3 + resi 234-266
select moving_domain_3, moving_domain_3 + resi 277-278
color pink, moving_domain_3

# Color bending residues
select bending_residues_1, resi 306-310
color green, bending_residues_1
select bending_residues_2, resi 267-273
color green, bending_residues_2
select bending_residues_3, resi 275-276
color green, bending_residues_3
select bending_residues_4, resi 279-280
color green, bending_residues_4

set dash_gap, 0
set dash_radius, 0.2

# === SCREW AXIS ARROWS ===
load 1yy9_a_epidermalgrowthfactor_A_1ivo_a_epidermalgrowthfactor_A_arrows.pdb
hide everything, 1yy9_a_epidermalgrowthfactor_A_1ivo_a_epidermalgrowthfactor_A_arrows

# Arrow 1: Domain 0 movement (116.3°)
select arrow_shaft_1, chain A and resn SHF and resi 100
select arrow_head_1, chain A and resn ARH and resi 120
show sticks, arrow_shaft_1
show sticks, arrow_head_1
color blue, arrow_shaft_1
color red, arrow_head_1
set stick_radius, 0.3, arrow_shaft_1
set stick_radius, 0.25, arrow_head_1
bond arrow_shaft_1, arrow_shaft_1
bond arrow_head_1, arrow_head_1

# Arrow 2: Domain 1 movement (30.8°)
select arrow_shaft_2, chain B and resn SHF and resi 150
select arrow_head_2, chain B and resn ARH and resi 170
show sticks, arrow_shaft_2
show sticks, arrow_head_2
color blue, arrow_shaft_2
color yellow, arrow_head_2
set stick_radius, 0.3, arrow_shaft_2
set stick_radius, 0.25, arrow_head_2
bond arrow_shaft_2, arrow_shaft_2
bond arrow_head_2, arrow_head_2

# Arrow 3: Domain 3 movement (32.9°)
select arrow_shaft_3, chain C and resn SHF and resi 200
select arrow_head_3, chain C and resn ARH and resi 220
show sticks, arrow_shaft_3
show sticks, arrow_head_3
color blue, arrow_shaft_3
color pink, arrow_head_3
set stick_radius, 0.3, arrow_shaft_3
set stick_radius, 0.25, arrow_head_3
bond arrow_shaft_3, arrow_shaft_3
bond arrow_head_3, arrow_head_3

set auto_bond, 0
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
print 'Moving domain 0: red, rotation 116.3°'
print 'Moving domain 1: yellow, rotation 30.8°'
print 'Moving domain 3: pink, rotation 32.9°'
