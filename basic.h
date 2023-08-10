# ifndef BASIC_DEF
# define BASIC_DEF

class Vector2
{
public:
    int u = 0;
    int v = 0;
};

class Vector3
{
public:
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vector3();
    Vector3(double a, double b, double c);
    Vector3 operator+(Vector3 b);
    Vector3 operator-(Vector3 b);
    Vector3 operator*(double b);
    bool operator==(Vector3 b);
    bool operator!=(Vector3 b);
    double dotProd(Vector3 b);
    double lenVec();
    Vector3 crossProd(Vector3 b);
};

Vector3 normalize(Vector3 a);

double getAngleCos(Vector3 a, Vector3 b);

class Color
{
public:
    double r = 0.0;
    double g = 0.0;
    double b = 0.0;

    Color();
    Color(double a, double y, double c);
    bool operator==(Color color);
    Color operator=(Color color);
    Color operator*(double a);
    Color operator+(Color a);
};

class Plane
{
public:
    double A;
    double B;
    double C;
    double D;
};

# endif