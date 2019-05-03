To run :
Extensions require openCV, hence it is needed to run.

1. Ensure openCV’s include in skeleton.cpp is set to the correct folder, current
   version should be the default location for openCV 3.4.1_5 on MacOS. For lab
   machines, this should just be ‘#include <opencv/cv.h>’.

2. Set ‘OPENCV_FLAGS’ to the correct flags in the makefile. For lab machines,
   this should be ‘/usr/lib64/libopencv_core.so.2.4 /usr/lib64/libopencv_highgui.so.2.4’.

3. From the terminal, make and run as usual.

4. To change the texture of the scene, in TestModelH.h, change the values of
   ‘int setting’ for walls and ‘int settingBoxes’ for boxes.
   0 = original colours, 1 = marble, 2 = metal grill, 3 = woven wood.
