# RayTracer-Source

### The code renders a 3D scene with 2 light sources and 3 transparent spheres within a striped box. 
![image](MyImage.png)
### How to run the code:
1. Use the Makefile to compile
2. The code takes in the name of input file(without suffix) as a cmdline argument. 
   Test command goes like:
   ./raytracer testFilePrefix(string) maxDepthOfReflection(int)
   e.g. ./raytracer  test_samp  2
3. The generated PPM file called MyImage by default.
    Contains a test sample: test_samp.txt
