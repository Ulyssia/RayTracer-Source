# include"basic.h"
# include<cmath>

Vector3::Vector3()
{
    x = x;
    y = y;
    z = z;
}

Vector3::Vector3(double a, double b, double c)
{
    x = a;
    y = b;
    z = c;
}

Vector3 Vector3::operator+(Vector3 b)
{
    Vector3 c;
    c.x = x+b.x;
    c.y = y+b.y;
    c.z = z+b.z;
    return c;
}

Vector3 Vector3::operator-(Vector3 b)
{
    Vector3 c;
    c.x = x-b.x;
    c.y = y-b.y;
    c.z = z-b.z;
    return c;
}

Vector3 Vector3::operator*(double b)
{
    Vector3 c{0,0,0};
    c.x = x*b;
    c.y = y*b;
    c.z = z*b;
    return c;
}

bool Vector3::operator==(Vector3 b)
{
    if(x == b.x&&y == b.y&&z == b.z)
    {
        return 1;
    }
    return 0;
}

bool Vector3::operator!=(Vector3 b)
{
    if(x != b.x||y != b.y||z != b.z)
    {
        return 1;
    }
    return 0;
}

double Vector3::dotProd(Vector3 b)
{
    return(x*b.x+y*b.y+z*b.z);
}

double Vector3::lenVec()
{
    return(sqrt(x*x+y*y+z*z));
}

Vector3 Vector3::crossProd(Vector3 b)
{
    Vector3 c{0,0,0};
    c.x = y*b.z-z*b.y;
    c.y = z*b.x-x*b.z;
    c.z = x*b.y-y*b.x;
    return c;
}

Vector3 normalize(Vector3 a)
{
    return(a*(1.0/a.lenVec()));
}

double getAngleCos(Vector3 a, Vector3 b)
{
    double len_a = a.lenVec();
    double len_b = b.lenVec();
    double temp = a.dotProd(b);
    double res = temp/(len_a*len_b);
    return res;
}

Color::Color()
{
    r=r;
    g=g;
    b=b;
}

Color::Color(double a, double y, double c)
{
    r = a;
    g = y;
    b = c;
}

bool Color::operator==(Color color)
{
    bool res = 0;
    if(r==color.r&&g==color.g&&b==color.b)
    {
        res = 1;
    }
    return res;
}

Color Color::operator=(Color color)
{
    Color res;
    res.r = color.r;
    res.g = color.g;
    res.b = color.b;
    return res;
}

Color Color::operator*(double a)
{
    Color res;
    res.r *= a;
    res.g *= a;
    res.b *= a;
    return res;
}

Color Color::operator+(Color a)
{
    Color res;
    res.r += a.r;
    res.g += a.g;
    res.b += a.b;
    return res;
}

long double min(long double a, long double b)
{
    if(a>b)
        return b;
    else
        return a;
}