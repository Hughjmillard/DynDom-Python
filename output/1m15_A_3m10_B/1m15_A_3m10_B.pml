reinitialize
load 1m15_A_3m10_B.pdb
bg_color white
color grey
select region0, resi 2-95
select region0, region0 + resi 124-152
select region0, region0 + resi 162-168
select region0, region0 + resi 190-206
select region0, region0 + resi 217-227
select region0, region0 + resi 256-274
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 97-124
select region1, region1 + resi 155-159
select region1, region1 + resi 160-162
select region1, region1 + resi 211-212
select region1, region1 + resi 230-250
select region1, region1 + resi 254-256
select region1, region1 + resi 276-311
select region1, region1 + resi 321-354
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 171-180
select region2, region2 + resi 210-211
select region2, region2 + resi 212-214
select region2, region2 + resi 215-217
set_color colour2 = [255,255,0  ]
color colour2, region2
select region3, resi 95-97
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 152-155
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 161-162
set_color colour5 = [0  ,255,0  ]
color colour5, region5
select region6, resi 227-230
set_color colour6 = [0  ,255,0  ]
color colour6, region6
select region7, resi 252-256
set_color colour7 = [0  ,255,0  ]
color colour7, region7
select region8, resi 274-276
set_color colour8 = [0  ,255,0  ]
color colour8, region8
select region9, resi 168-171
set_color colour9 = [0  ,255,0  ]
color colour9, region9
select region10, resi 182-190
set_color colour10 = [0  ,255,0  ]
color colour10, region10
select region11, resi 206-210
set_color colour11 = [0  ,255,0  ]
color colour11, region11
select region12, resi 216-217
set_color colour12 = [0  ,255,0  ]
color colour12, region12
set dash_gap, 0
set dash_radius, 0.2
