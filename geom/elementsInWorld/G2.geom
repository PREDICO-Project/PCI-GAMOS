############################
##   G2-grating phantom   ##
############################


#Params definition
:P G2Position 1900.
:P G2Period 0.006
:P G2DutyCycle 0.5
:P G2Thickness 0.1
:P G2SizeX 50.
:P G2SizeY 50.

:P NumberCopiesG2 $G2SizeX/$G2Period 

########################
##   VOL definition   ##
########################

#layer of the grating made by Silicon
:VOLU G2 BOX 0.5*$G2SizeX 0.5*$G2SizeY 0.5*$G2Thickness G4_Si
:PLACE G2 1 world RM2 $G2POSX 0. $G2Position

# Gold part 
:VOLU G2_Au BOX 0.5*$G2Period*$G2DutyCycle 0.5*$G2SizeY 0.5*$G2Thickness G4_Au
:COLOR G2_Au 1.0 0. 0. 1.

:PLACE_PARAM G2_Au 2 G2 LINEAR_X RM0 $NumberCopiesG2 $G2Period -0.5*$G2SizeX+0.5*$G2Period*$G2DutyCycle
