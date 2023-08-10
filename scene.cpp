# include "basic.h"
# include "scene.h"

double getAttF(Light light, double d)
{
    return(1.0/(light.attParam.x + light.attParam.y*d + light.attParam.z*d*d));
}

bool Sphere::operator==(Sphere s)
{
    if(center == s.center&&radius == s.radius)
    {
        return 1;
    }
    return 0;
}

bool Triangle_Face::operator==(Triangle_Face t)
{
    if(t.p1==p1&&t.p2==p2&&t.p3==p3)
    {
        return 1;
    }
    return 0;
}

Vector3 getFaceNormal(Triangle_Face tf, Vector3* vertex)
{
    Vector3 temp1 = vertex[tf.p2] - vertex[tf.p1];
    Vector3 temp2 = vertex[tf.p3] - vertex[tf.p1];
    Vector3 norm = temp1.crossProd(temp2);
    norm = normalize(norm);
    return norm;
}