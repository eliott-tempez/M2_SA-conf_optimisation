load ./saconf_out/PDB/3QKD.pdb, 3QKD 
color grey, 3QKD
select protch, chain a, 3QKD
color cyan, protch
hide line
show cartoon, 3QKD
show surface, protch
remove resname hoh
set transparency, 0.3
select ligand, het 
show sticks, ligand
util.cbay("ligand")
select mut, chain a and resid 24+25+26+83+84+92+97+109+113+137+142+146+158+168+189 
show sticks, mut
select consPos, chain a and  resid 137+187 
color magenta, consPos
select varPos_NoSS, chain a and  resid 16+22+96+116+118+145+147+148+169+170+171+172+173+176+188 
color slate, varPos_NoSS
select varPos_withSS, chain a and  resid 4+5+6+7+8+9+10+11+12+13+14+15+17+18+19+20+23+24+25+26+83+84+85+86+87+88+89+90+91+92+93+94+95+97+98+99+100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+117+119+120+121+122+123+124+125+126+127+128+129+130+131+133+135+138+139+140+141+142+143+144+149+151+152+153+154+155+156+157+159+160+161+162+163+164+165+166+167+168+174+175+180+181+182+183+184+185+190+191+192+193+194+195 
color blue, varPos_withSS
bg_colour white 
set surface_quality, 2 
set cartoon_fancy_helices, 1 
set cartoon_fancy_sheets, 1 
set stick_radius, 0.2 
