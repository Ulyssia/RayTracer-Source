# ifndef SCENE_DEF
# define SCENE_DEF
# include "basic.h"

class MtlColor:public Color
{
public:
    double Osr = 0.0;
    double Osg = 0.0;
    double Osb = 0.0;
    double ka = 0.0;
    double kd = 0.0;
    double ks = 0.0;
    double exp_n = 0.0;
    double alpha = 0.0;
    double eta = 0.0;
};

class Light
{
public:
    Vector3 pos;
    double w = 2.0;
    Color light_color = {1.0,1.0,1.0};
    bool isAtt = 0;
    Vector3 attParam;  
};

double getAttF(Light light, double d);

class Camera                                                                                                                               
{
public:
    Vector3 eye;
    Vector3 viewdir;
    Vector3 updir;
    double vfov;
};

class Ray
{
public:
    Vector3 start_point;
    Vector3 dir;
    Color color = {0,0,0};
    bool isInside = 0;
};

class Sphere
{
public:
    Vector3 center;
    double radius;
    MtlColor surf_material;
    bool useText = false;

    bool operator==(Sphere s);
};

class Triangle_Face
{
public:
    int p1 = 0;
    int p2 = 0;
    int p3 = 0;
    bool isFlat = true;
    int normal1 = 0;
    int normal2 = 0;
    int normal3 = 0;
    bool useText = false;
    int vt1 = 0;
    int vt2 = 0;
    int vt3 = 0;
    MtlColor material;

    bool operator==(Triangle_Face t);
};

Vector3 getFaceNormal(Triangle_Face tf, Vector3* vertex);

class Scene_objs
{
public:
    Sphere* sphere;
    Triangle_Face* triangle;
    Vector3* vertex;
    Vector3* v_normal;
    Vector2* v_coord;
    Vector2 text_size;
    Color** texture;
    Color background;
};

# endif