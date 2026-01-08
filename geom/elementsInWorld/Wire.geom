//*******VOLUMES DEFINITION******//
//Rotation Matrix
:ROTM RM90 -90. 0. 0.

//*******PARAMS**********//

//Test Object
:P InnerRadius 2
//:P OuterRadius 0.0125
:P OuterRadius 5
:P HalfLengthZ 40
//********TEST OBJECT**********//

//Solid water sphere filled with a teflon box
//:VOLU WIRE TUBE $InnerRadius $OuterRadius $HalfLengthZ "AAPM_PMMA"
:VOLU WIRE TUBE $InnerRadius $OuterRadius $HalfLengthZ G4_WATER
:PLACE WIRE 1 world RM90 0. 0. $SOS
:COLOUR WIRE 1.0 0. 0.
