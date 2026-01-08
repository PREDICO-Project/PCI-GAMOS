######################## 
##    WORLD SETUP     ## 
######################## 

#Materials 
:include geom/MyMaterials.txt 

#Rot Matrix

:ROTM RM0 0. 0. 0.
:ROTM RM1 0. 0. 0. 
:ROTM RM2 0. 0. 0.

#Params definition 
#:P SOS 1200.0

# Source-G1 distance 
:P SG1D 403.2
:P SDD 2000.0

:P WorldDimen $SDD*2+50 

# G2/G2 shift
:P G1POSX 0.
:P G2POSX 0.

#World shape and composition 
:VOLU world BOX 0.5*$WorldDimen 0.5*$WorldDimen 0.5*$WorldDimen G4_AIR 
#:VIS world OFF

#Jaws 
#:include geom/jaws.geom 

#Phantom (If we use a DICOM object the following line must be commented or empty) 

#Detector 
#:include geom/elementsInWorld/Cylinder.geom
#:include geom/Wire.geom
#:include geom/Wedge.geom
#:include geom/elementsInWorld/G1.geom
#:include /elementsInWorld/geom/G2.geom
#:include geom/detector.geom
:include geom/elementsInWorld/Grating_test.geom
#:include geom/elementsInWorld/Abs_Grating.geom
