reinitialize
load 1adg_A_2ohx_B.pdb
bg_color white
color grey
select region0, resi 4-174
select region0, region0 + resi 290-291
select region0, region0 + resi 319-372
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 177-282
select region1, region1 + resi 286-287
select region1, region1 + resi 297-311
select region1, region1 + resi 315-319
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 174-177
set_color colour2 = [0  ,255,0  ]
color colour2, region2
select region3, resi 286-290
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 291-297
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 315-319
set_color colour5 = [0  ,255,0  ]
color colour5, region5
set dash_gap, 0
set dash_radius, 0.2
