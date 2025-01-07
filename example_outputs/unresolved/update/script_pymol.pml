load example_outputs/unresolved/update/PDB/3QYM.pdb, 3QYM 
color grey, 3QYM
select protch, chain b, 3QYM
color cyan, protch
hide line
show cartoon, 3QYM
show surface, protch
remove resname hoh
set transparency, 0.3
select ligand, het 
show sticks, ligand
util.cbay("ligand")
select mut, chain b and resid 144 
show sticks, mut
select consPos, chain b and  resid 128+129+130+131+134+137+138+139+140+141+144+153+154+155+160+161+164+167+168+169+170+173+175+176+177+178+179+180+182+187+188+190+191+192+193+194+195+199+200+201+202+203+205+206+214+215+220+221+224+225+227+228+229+230+231+232+233+235+237+238+239+240+242+244+246+247+248+249+250+251+252+259+261+262+265+267+268+269+271+273+275+279+280+281+284+285+289+291+293+294+295+296+297+298+299+300+301+303+307+308+309+311+312 
color magenta, consPos
select varPos_NoSS, chain b and  resid 133+135+142+143+145+150+152+156+157+162+163+166+172+174+181+183+185+189+196+197+198+204+210+211+212+213+216+217+218+219+222+223+226+234+236+241+255+256+257+260+263+264+266+276+277+278+282+283+286+288+315+316+317+318 
color slate, varPos_NoSS
select varPos_withSS, chain b and  resid 132+136+146+151+158+165+207+208+209+253+270+272+274+290+292+302+304+305+306+313+314 
color blue, varPos_withSS
bg_colour white 
set surface_quality, 2 
set cartoon_fancy_helices, 1 
set cartoon_fancy_sheets, 1 
set stick_radius, 0.2 
