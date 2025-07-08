reinitialize
load 1adg_A_2ohx_B.pdb
bg_color white
color grey
select region0, resi 1-171
select region0, region0 + resi 287-288
select region0, region0 + resi 316-369
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 174-282
select region1, region1 + resi 286-287
select region1, region1 + resi 294-311
select region1, region1 + resi 315-316
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 171-174
set_color colour2 = [0  ,255,0  ]
color colour2, region2
select region3, resi 283-287
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 288-294
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 312-316
set_color colour5 = [0  ,255,0  ]
color colour5, region5
set dash_gap, 0
set dash_radius, 0.2
