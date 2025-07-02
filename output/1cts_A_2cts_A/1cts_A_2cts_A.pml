reinitialize
load 1cts_A_2cts_A.pdb
bg_color white
color grey
select region0, resi 1-53
select region0, region0 + resi 65-269
select region0, region0 + resi 331-336
select region0, region0 + resi 377-434
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 55-62
select region1, region1 + resi 273-276
select region1, region1 + resi 280-331
select region1, region1 + resi 338-340
select region1, region1 + resi 345-372
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 53-55
set_color colour2 = [0  ,255,0  ]
color colour2, region2
select region3, resi 63-65
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 269-273
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 277-278
set_color colour5 = [0  ,255,0  ]
color colour5, region5
select region6, resi 279-280
set_color colour6 = [0  ,255,0  ]
color colour6, region6
select region7, resi 336-338
set_color colour7 = [0  ,255,0  ]
color colour7, region7
select region8, resi 341-345
set_color colour8 = [0  ,255,0  ]
color colour8, region8
select region9, resi 373-375
set_color colour9 = [0  ,255,0  ]
color colour9, region9
select region10, resi 376-377
set_color colour10 = [0  ,255,0  ]
color colour10, region10
set dash_gap, 0
set dash_radius, 0.2
