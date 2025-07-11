reinitialize
load 1yy9_a_epidermalgrowthfactor_A_1ivo_a_epidermalgrowthfactor_A.pdb
bg_color white
color grey
select region0, resi 274-274
select region0, region0 + resi 281-305
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 311-508
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 2-5
select region2, region2 + resi 7-228
select region2, region2 + resi 230-230
select region2, region2 + resi 232-233
set_color colour2 = [255,255,0  ]
color colour2, region2
select region3, resi 6-6
select region3, region3 + resi 229-229
select region3, region3 + resi 231-231
select region3, region3 + resi 234-266
select region3, region3 + resi 277-278
set_color colour3 = [255,100,255]
color colour3, region3
select region4, resi 306-310
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 267-273
set_color colour5 = [0  ,255,0  ]
color colour5, region5
select region6, resi 275-276
set_color colour6 = [0  ,255,0  ]
color colour6, region6
select region7, resi 279-280
set_color colour7 = [0  ,255,0  ]
color colour7, region7
set dash_gap, 0
set dash_radius, 0.2
