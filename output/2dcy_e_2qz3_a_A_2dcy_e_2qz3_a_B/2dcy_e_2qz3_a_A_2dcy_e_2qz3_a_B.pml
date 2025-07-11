reinitialize
load 2dcy_e_2qz3_a_A_2dcy_e_2qz3_a_B.pdb
bg_color white
color grey
select region0, resi 2-40
select region0, region0 + resi 48-64
select region0, region0 + resi 80-91
select region0, region0 + resi 129-141
select region0, region0 + resi 150-156
select region0, region0 + resi 166-181
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 68-77
select region1, region1 + resi 94-97
select region1, region1 + resi 102-126
select region1, region1 + resi 144-147
select region1, region1 + resi 159-163
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 41-47
set_color colour2 = [0  ,255,0  ]
color colour2, region2
select region3, resi 65-67
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 78-79
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 92-93
set_color colour5 = [0  ,255,0  ]
color colour5, region5
select region6, resi 98-101
set_color colour6 = [0  ,255,0  ]
color colour6, region6
select region7, resi 127-128
set_color colour7 = [0  ,255,0  ]
color colour7, region7
select region8, resi 142-143
set_color colour8 = [0  ,255,0  ]
color colour8, region8
select region9, resi 148-149
set_color colour9 = [0  ,255,0  ]
color colour9, region9
select region10, resi 157-158
set_color colour10 = [0  ,255,0  ]
color colour10, region10
select region11, resi 164-165
set_color colour11 = [0  ,255,0  ]
color colour11, region11
set dash_gap, 0
set dash_radius, 0.2
