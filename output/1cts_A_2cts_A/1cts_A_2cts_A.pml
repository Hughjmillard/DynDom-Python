reinitialize
load 1cts_A_2cts_A.pdb
bg_color white
color grey
select region0, resi 3-55
select region0, region0 + resi 67-271
select region0, region0 + resi 333-338
select region0, region0 + resi 379-436
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 57-62
select region1, region1 + resi 64-65
select region1, region1 + resi 275-276
select region1, region1 + resi 277-278
select region1, region1 + resi 280-281
select region1, region1 + resi 282-333
select region1, region1 + resi 347-372
select region1, region1 + resi 374-375
select region1, region1 + resi 377-378
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 55-57
set_color colour2 = [0  ,255,0  ]
color colour2, region2
select region3, resi 65-67
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 271-275
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 279-280
set_color colour5 = [0  ,255,0  ]
color colour5, region5
select region6, resi 281-282
set_color colour6 = [0  ,255,0  ]
color colour6, region6
select region7, resi 338-340
set_color colour7 = [0  ,255,0  ]
color colour7, region7
select region8, resi 343-347
set_color colour8 = [0  ,255,0  ]
color colour8, region8
select region9, resi 375-377
set_color colour9 = [0  ,255,0  ]
color colour9, region9
select region10, resi 378-379
set_color colour10 = [0  ,255,0  ]
color colour10, region10
set dash_gap, 0
set dash_radius, 0.2
