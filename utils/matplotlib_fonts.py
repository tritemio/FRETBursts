# encoding: utf-8  μ αβγδ εζηθ κλμν ξοπρ ςστυ φχψω

##
# Experiments with fonts
#
# NOTE: to change default font in realtime change font.family to your
# font name. Keeping the family to sans-serif and changing font.sans-serif
# change the fornt only the first time!!!
#

font1 = {'fontname':'Liberation Sans','fontsize':16} # best general font
font2 = {'fontname':'Arial','fontsize':16} # Best for numbers, worst greek
font3 = {'fontname':'Dejavu Sans','fontsize':16}
font4 = {'fontname':'Ubuntu','fontsize':16} # Nice 1 number very sim. to Arial
font5 = {'fontname':'cmex10', 'fontsize':16}
font5 = {'fontname':'Computer Modern Sans serif', 'fontsize':16}

rcParams["font.family"] = font5['fontname']

greek = u"ταβγδεζηθικλμνξοπρςστυφχψω"
greek_math_it = u"𝛼𝛽𝛾𝛿𝜀𝜁𝜂𝜃𝜆𝜇𝜈𝜉𝜋𝜌𝜍𝜎𝜏𝜐𝜑𝜒𝜓𝜔𝜕𝜖𝜗𝜘𝜙𝜚" # apparently does not work
s = " 0123456789 "

figure()
title(s) # -> uses default font

# Use custom fonts
text(0.1,0.9, font1['fontname']+s+greek, **font1) 
text(0.1,0.7, font2['fontname']+s+greek, **font2)
text(0.1,0.3, font3['fontname']+s+greek, **font3)
text(0.1,0.1, font4['fontname']+s+greek, **font4)
text(0.1,0.5, font5['fontname']+s+greek, **font4)

## Last (default) setting before show() takes effect
#rcParams["font.family"] = font1['fontname']

show()

