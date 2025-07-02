reinitialize
load calmodulin1_A_calmodulin2_B.pdb
bg_color white
color grey
select region0, resi 5-68
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 74-75
select region1, region1 + resi 76-143
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 68-76
set_color colour2 = [0  ,255,0  ]
color colour2, region2
set dash_gap, 0
set dash_radius, 0.2
