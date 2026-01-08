########################
##  CYLINDER phantom  ##
########################

##Rotation Matrix
:ROTM RM90 -90. 0. 0.
#Params definition
#:P InnerRadius 0.*mm
#:P HalfLengthZ 0.8*mm
#:P OuterRadius 0.100*mm

:P InnerRadius 1*mm
:P OuterRadius 2*mm

#:P InnerRadius 0*mm
#:P OuterRadius 1*mm
:P HalfLengthZ 8*mm

########################
##   VOL definition   ##
########################

#Solid PMMA cylinder 
:VOLU CYL TUBE $InnerRadius $OuterRadius $HalfLengthZ G4_WATER
:PLACE CYL 1 world RM90 0. 0. $SOD
