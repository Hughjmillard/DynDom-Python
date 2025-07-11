reinitialize
load 1m15a_A_3m10b_B.pdb
bg_color white
color grey
select region0, resi 2-94
select region0, region0 + resi 125-151
select region0, region0 + resi 163-167
select region0, region0 + resi 190-205
select region0, region0 + resi 218-226
select region0, region0 + resi 257-273
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 98-122
select region1, region1 + resi 156-160
select region1, region1 + resi 211-211
select region1, region1 + resi 231-251
select region1, region1 + resi 277-310
select region1, region1 + resi 321-353
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 172-181
select region2, region2 + resi 212-215
set_color colour2 = [255,255,0  ]
color colour2, region2
select region3, resi 95-97
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 123-124
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 152-155
set_color colour5 = [0  ,255,0  ]
color colour5, region5
select region6, resi 161-162
set_color colour6 = [0  ,255,0  ]
color colour6, region6
select region7, resi 227-230
set_color colour7 = [0  ,255,0  ]
color colour7, region7
select region8, resi 252-256
set_color colour8 = [0  ,255,0  ]
color colour8, region8
select region9, resi 274-276
set_color colour9 = [0  ,255,0  ]
color colour9, region9
select region10, resi 168-171
set_color colour10 = [0  ,255,0  ]
color colour10, region10
select region11, resi 182-189
set_color colour11 = [0  ,255,0  ]
color colour11, region11
select region12, resi 206-210
set_color colour12 = [0  ,255,0  ]
color colour12, region12
select region13, resi 216-217
set_color colour13 = [0  ,255,0  ]
color colour13, region13
set dash_gap, 0
set dash_radius, 0.2
