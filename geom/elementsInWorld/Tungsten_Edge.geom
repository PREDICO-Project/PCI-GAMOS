//*******VOLUMES DEFINITION******//
//Rotation Matrix
:ROTM RM02z 0. 0. 4.5

//*******PARAMS**********//

//Test Object

:P EdgeSideX 75.
:P EdgeSideY 75.

:P EdgeWidth 1.


//********TEST OBJECT**********//

//Tungsten Edge for MTF calculation
:VOLU EDGE BOX 0.5*$EdgeSideX 0.5*$EdgeSideY 0.5*$EdgeWidth G4_W

:PLACE EDGE 1 world RM02z 0. 0. $SOD-0.5*$EdgeWidth

:COLOUR EDGE 1.0 0. 0.
