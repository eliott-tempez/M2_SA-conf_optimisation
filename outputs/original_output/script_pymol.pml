load ./saconf_out/PDB/1HHP.pdb, 1HHP 
color grey, 1HHP
select protch, chain a, 1HHP
color cyan, protch
hide line
show cartoon, 1HHP
show surface, protch
remove resname hoh
set transparency, 0.3
select ligand, het 
show sticks, ligand
util.cbay("ligand")
select mut, chain a and resid 3+4+7+9+10+11+12+13+14+15+16+17+19+20+22+31+32+33+34+35+36+37+39+40+41+42+43+46+47+55+56+57+58+60+61+62+63+66+67+68+69+71+72+73+75+76+77+79+82+85+89+92+93+95+96+99 
show sticks, mut
select consPos, chain a and  resid 7+28+30+68+82+84 
color magenta, consPos
select varPos_NoSS, chain a and  resid 6+18+32+35+40+41+43+46+50+51+64+67+69+70+75+77+79+87+89+96 
color slate, varPos_NoSS
select varPos_withSS, chain a and  resid 3+4+5+10+13+14+15+16+21+22+23+24+25+31+34+36+37+38+39+42+44+47+48+49+54+57+58+59+60+62+63+66+72+73+74+78+80+83+86+90+91+92+93+95+97+98 
color blue, varPos_withSS
bg_colour white 
set surface_quality, 2 
set cartoon_fancy_helices, 1 
set cartoon_fancy_sheets, 1 
set stick_radius, 0.2 
