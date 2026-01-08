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
#:P SOD 400.0*mm
:P SDD 2000.0*mm 
:P WorldDimen $SDD*2+50 
:P SOD 400.8*mm

:P G1POSX 0.

#World shape and composition 
:VOLU world BOX 0.5*$WorldDimen 0.5*$WorldDimen 0.5*$WorldDimen G4_Galactic

#Detector 
#:include geom/elementsInWorld/Cylinder.geom
:include geom/elementsInWorld/Sphere.geom
#:include geom/elementsInWorld/test_detector.geom
#:include geom/elementsInWorld/Grating_test.geom
#:include geom/elementsInWorld/double_slit.geom
#:include geom/elementsInWorld/Abs_Grating.geom
