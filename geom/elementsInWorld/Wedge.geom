//*******VOLUMES DEFINITION******//
//Rotation Matrix
:ROTM RMTrap 90. 0. 0.

//*******PARAMS**********//

//Test Object
:P dz 5 // Depth (height in Z)
:P dy1 5 // Half height of the bottom face
:P dx1 10 // Half width of the bottom face
:P dx2 0.0000001 // Half width of the top face
:P alpha1 0. // Bottom face tilt angle
:P dy2 5 // Half height of the top face
:P dx3 10. // Half width of the top face (upper plane)
:P dx4 0.0000001 // Half width of the bottom face (upper plane) 
:P alpha2 0.0 // Top face tilt plane
//********TEST OBJECT**********//

//Solid water sphere filled with a teflon box

:VOLU WEDGE TRAP $dz $alpha1 $alpha2 $dy1 $dx1 $dx2 0.0 $dy2 $dx3 $dx4 0.0 G4_WATER
:PLACE WEDGE 1 world RMTrap 0. 0. $SOS
:COLOUR WEDGE 1.0 0. 0.


//:VOLU NAME TRAP 
//•Half-length along the z-axis (=pDz)
//• Polar angle of the line joining the centres of the faces at -/+pDz
//• Azimuthal angle of the line joining the centre of the face at -pDzto the centre of the face at +pDz
//• Half-length along y of the face at -pDz (=pDy1)
//• Half-length along x of the side at y=-pDy1 of the face at -pDz
//• Half-length along x of the side at y=+pDy1 of the face at -pDz
//• Angle with respect to the y axis from the centre of the side at y=-pDy1 to the centre at y=+pDy1 of the face at -pDz
//• Half-length along y of the face at +pDz (=pDy2)
//• Half-length along x of the side at y=-pDy2 of the face at +pDz
//• Half-length along x of the side at y=+pDy2 of the face at +pDz
//• Angle with respect to the y axis from the centre of the side at y=-pDy2 to the centre at y=+pDy2 of the face at +pDz
