reinitialize
load 4ake_A_2eck_B.pdb
bg_color white
color grey
select region0, resi 1-25
select region0, region0 + resi 63-111
select region0, region0 + resi 169-210
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 116-152
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 29-58
set_color colour2 = [255,255,0  ]
color colour2, region2
select region3, resi 112-115
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 153-168
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 26-28
set_color colour5 = [0  ,255,0  ]
color colour5, region5
select region6, resi 59-62
set_color colour6 = [0  ,255,0  ]
color colour6, region6
set dash_gap, 0
set dash_radius, 0.2
