reinitialize
load thior1r5m1_A_thior1r5m2_B.pdb
bg_color white
color grey
select region0, resi 1-108
select region0, region0 + resi 244-313
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 113-238
select region1, region1 + resi 243-244
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 108-113
set_color colour2 = [0  ,255,0  ]
color colour2, region2
select region3, resi 239-244
set_color colour3 = [0  ,255,0  ]
color colour3, region3
set dash_gap, 0
set dash_radius, 0.2
