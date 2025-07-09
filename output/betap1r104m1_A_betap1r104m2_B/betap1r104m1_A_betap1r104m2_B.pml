reinitialize
load betap1r104m1_A_betap1r104m2_B.pdb
bg_color white
color grey
select region0, resi 78-88
select region0, region0 + resi 112-128
select region0, region0 + resi 135-137
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 14-74
select region1, region1 + resi 77-78
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 1-6
select region2, region2 + resi 90-108
select region2, region2 + resi 111-112
select region2, region2 + resi 130-132
select region2, region2 + resi 134-135
select region2, region2 + resi 140-216
set_color colour2 = [255,255,0  ]
color colour2, region2
select region3, resi 11-14
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 75-78
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 7-14
set_color colour5 = [0  ,255,0  ]
color colour5, region5
select region6, resi 88-90
set_color colour6 = [0  ,255,0  ]
color colour6, region6
select region7, resi 109-112
set_color colour7 = [0  ,255,0  ]
color colour7, region7
select region8, resi 128-130
set_color colour8 = [0  ,255,0  ]
color colour8, region8
select region9, resi 133-135
set_color colour9 = [0  ,255,0  ]
color colour9, region9
select region10, resi 137-140
set_color colour10 = [0  ,255,0  ]
color colour10, region10
set dash_gap, 0
set dash_radius, 0.2
