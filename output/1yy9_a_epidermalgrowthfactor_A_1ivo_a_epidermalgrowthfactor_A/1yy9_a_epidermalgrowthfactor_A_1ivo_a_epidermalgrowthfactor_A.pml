reinitialize
load 1yy9_a_epidermalgrowthfactor_A_1ivo_a_epidermalgrowthfactor_A.pdb
bg_color white
color grey
select region0, resi 273-276
select region0, region0 + resi 280-306
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 310-509
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 2-6
select region2, region2 + resi 7-229
select region2, region2 + resi 230-231
select region2, region2 + resi 232-234
set_color colour2 = [255,255,0  ]
color colour2, region2
select region3, resi 6-7
select region3, region3 + resi 229-230
select region3, region3 + resi 231-232
select region3, region3 + resi 234-265
select region3, region3 + resi 271-273
select region3, region3 + resi 276-280
set_color colour3 = [255,100,255]
color colour3, region3
select region4, resi 306-310
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 267-273
set_color colour5 = [0  ,255,0  ]
color colour5, region5
set dash_gap, 0
set dash_radius, 0.2
