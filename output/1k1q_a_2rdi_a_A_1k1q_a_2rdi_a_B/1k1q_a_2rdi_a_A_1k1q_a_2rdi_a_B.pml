reinitialize
load 1k1q_a_2rdi_a_A_1k1q_a_2rdi_a_B.pdb
bg_color white
color grey
select region0, resi 1-34
select region0, region0 + resi 46-117
select region0, region0 + resi 118-235
select region0, region0 + resi 238-239
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 239-279
select region1, region1 + resi 280-330
select region1, region1 + resi 331-336
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 235-238
set_color colour2 = [0  ,255,0  ]
color colour2, region2
set dash_gap, 0
set dash_radius, 0.2
