reinitialize
load 1cll_a_from_dyndomserver_A_1cll_a_screw1_A.pdb
bg_color white
color grey
select region0, resi 5-98
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 98-143
set_color colour1 = [255,0  ,0  ]
color colour1, region1
set dash_gap, 0
set dash_radius, 0.2
