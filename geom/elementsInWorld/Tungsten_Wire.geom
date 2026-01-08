//*******VOLUMES DEFINITION******//
//Rotation Matrix
:ROTM RM90 -90. -88. 0.

//*******PARAMS**********//

//Test Object
:P InnerRadius 0.
//:P OuterRadius 0.0125
:P OuterRadius 0.1
:P HalfLengthZ 40
//********TEST OBJECT**********//

//Solid water sphere filled with a teflon box
:VOLU WIRE TUBE $InnerRadius $OuterRadius $HalfLengthZ G4_W
:PLACE WIRE 1 world RM90 0. 0. $SOD
:COLOUR WIRE 1.0 0. 0.
