########################
##   SPHERE phantom   ##
########################

#Params definition
#:P SphereRadius 0.15*mm
:P SphereRadius 0.2*mm

########################
##   VOL definition   ##
########################

#Solid water sphere 
:VOLU test ORB $SphereRadius G4_WATER
:PLACE test 1 world RM1 0. 0. $SOD

