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
:P SDD 2000.0*mm
# Source-G1 distance
:P SOD 401*mm
:P SG1D 403.2*mm 

:P WorldDimen $SDD*2+50 

# G2 shift
:P G2POSX 0.
:P G1POSX 0.

#World shape and composition 
:VOLU world BOX 0.5*$WorldDimen 0.5*$WorldDimen 0.5*$WorldDimen G4_Galactic
#:VIS world OFF

#Detector 
:include geom/elementsInWorld/Cylinder.geom
:include geom/elementsInWorld/Grating_test.geom
#:include geom/elementsInWorld/Sphere.geom
#:include geom/elementsInWorld/Wire_fibers.geom
#:include geom/elementsInWorld/Wedge.geom
#:include geom/elementsInWorld/detector.geom
