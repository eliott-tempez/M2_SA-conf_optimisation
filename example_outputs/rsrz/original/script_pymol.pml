load ./saconf_out/PDB/1SDT.pdb, 1SDT 
color grey, 1SDT
select protch, chain a, 1SDT
color cyan, protch
hide line
show cartoon, 1SDT
show surface, protch
remove resname hoh
set transparency, 0.3
select ligand, het 
show sticks, ligand
util.cbay("ligand")
select mut, chain a and resid 3+4+10+13+14+16+19+20+22+31+32+33+34+35+36+37+39+40+41+42+43+46+47+55+56+57+58+60+61+62+63+66+67+68+69+71+72+73+75+76+77+79+82+85+89+92+93+95+96+99 
show sticks, mut
select consPos, chain a and  resid 6+7+8+9+11+12+17+19+20+21+24+25+26+27+28+29+30+33+45+52+53+55+56+57+61+62+65+68+71+76+77+79+80+81+82+84+85+88+94+97 
color magenta, consPos
select varPos_NoSS, chain a and  resid 3+5+10+13+14+15+16+18+22+23+32+35+40+41+43+46+47+48+50+51+54+58+59+60+64+66+67+69+70+72+74+75+86+87+89+91+92+93+95+96+98 
color slate, varPos_NoSS
select varPos_withSS, chain a and  resid 4+31+34+36+37+38+39+42+44+49+63+73+78+83+90 
color blue, varPos_withSS
bg_colour white 
set surface_quality, 2 
set cartoon_fancy_helices, 1 
set cartoon_fancy_sheets, 1 
set stick_radius, 0.2 
