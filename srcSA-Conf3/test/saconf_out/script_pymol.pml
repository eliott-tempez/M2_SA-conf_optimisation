load ./saconf_out/PDB/1HSI.pdb, 1HSI 
color grey, 1HSI
select protch, chain a, 1HSI
color cyan, protch
hide line
show cartoon, 1HSI
show surface, protch
remove resname hoh
set transparency, 0.3
select ligand, het 
show sticks, ligand
util.cbay("ligand")
select consPos, chain a and  resid 7+8+9+10+11+12+17+19+20+24+25+26+27+28+29+30+33+35+41+45+53+56+58+61+67+68+71+72+73+74+76+81+82+84+85+87+88+89+90+94+95+97 
color magenta, consPos
select varPos_NoSS, chain a and  resid 3+5+6+13+14+16+18+22+23+32+34+36+38+39+40+42+43+46+47+48+49+50+51+52+54+57+59+60+62+64+65+66+69+70+75+77+79+86+91+92+93+96+98 
color slate, varPos_NoSS
select varPos_withSS, chain a and  resid 4+15+21+31+37+44+55+63+78+80+83 
color blue, varPos_withSS
bg_colour white 
set surface_quality, 2 
set cartoon_fancy_helices, 1 
set cartoon_fancy_sheets, 1 
set stick_radius, 0.2 
