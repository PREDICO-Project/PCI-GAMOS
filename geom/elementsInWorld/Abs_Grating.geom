############################
##   G1-grating phantom   ##
############################



#Params definition

:P G1Position 403.2
:P G1Period 0.010
:P G1DutyCycle 0.5
:P G1Thickness 0.1
:P G1SizeX 50.
:P G1SizeY 50.

:P NumberCopiesG1 $G1SizeX/$G1Period 

########################
##   VOL definition   ##
########################

#layer of the grating made by Silicon
:VOLU G1 BOX 0.5*$G1SizeX 0.5*$G1SizeY 0.5*$G1Thickness G4_Pb
:PLACE G1 1 world RM1 $G1POSX 0. $G1Position

# Ni part 
:VOLU G1_Si BOX 0.5*$G1Period*(1-$G1DutyCycle) 0.5*$G1SizeY 0.5*$G1Thickness G4_Galactic
:COLOR G1_Si 1.0 0. 0. 1.

:PLACE_PARAM G1_Si 2 G1 LINEAR_X RM1 $NumberCopiesG1 $G1Period -0.5*$G1SizeX+0.5*$G1Period*$G1DutyCycle
