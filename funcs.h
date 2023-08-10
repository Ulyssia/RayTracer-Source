# ifndef FUNCS_DEF
# define FUNCS_DEF
# include "basic.h" 
# include "scene.h"

//get the specific point where the ray hit the object
Vector3 getPoint(Ray ray, double t);

//Decide the plane equation of a triangle
Plane getTriPlane(Triangle_Face tf, Vector3* vertex);

//Solve the intersection of a ray and a plane
double rayPlaneIntersect(Ray ray, Plane plane);

//Solve the barycentric coordinate of a given point relative to a triangle
Vector3 barycentric(Triangle_Face tf, Vector3* vertex, Vector3 p);

//Decide if ray/plane intersection is within triangle
bool isInTriangle(Vector3 p, Triangle_Face tf, Vector3* vertex);

//Decide the closest&positive ray/sphere intersection point
double closestIntersect(Sphere sphere_a, Ray ray);

//Decide color of a pixel given by float type position using Bilinear Interpolation
Color biLinInterplt(Scene_objs scene, double u, double v);

Color nesarestNeighbor(Scene_objs scene, double u, double v);

//Define the mapping from texture map to geometirc surfaces
Color textMap_sph(Scene_objs scene, Sphere sphere, Vector3 intersect);
Color textMap_Tri(Scene_objs scene, Triangle_Face tf, Vector3 intersect);

//get vector I out of incident ray for reflection
Ray getIncident(Ray init, Vector3 intersect);

Ray getReflect(Ray inc, Vector3 normal);

double getF0(MtlColor a);

double getFresnel(double f0, double cos_in);

Ray getTransmit(Ray init, Vector3 normal, double indx_i, double indx_t);

//Decide the color of pixels wrt illumination
Color ShadeRay(Light* lightSrc, Scene_objs scene, Camera camera, Vector3 view_init, Sphere sphere_a, Vector3 intersect);

Color ShadeRayTri(Light* lightSrc, Scene_objs scene, Camera camera, Triangle_Face tf, Vector3 intersect);

Color TraceSingle(Light* lightSrc, Scene_objs scene, Camera camera, int maxDepth, int depth, Ray cur_ray);

//Decide color for each pixel in the image
//For each pixel in the image, generate a ray that starts from the eye and goes across the pixel.
void RayTracer(Color** renderPic, int width, int height, Color bkg_color, Camera camera, Light* lightSrc, Scene_objs scene, int maxLevel, Vector3 ul, Vector3 dh, Vector3 dv);

void clamp(Color** renderPic, int width, int height);

# endif