# define _USE_MATH_DEFINES
# include "funcs.h"
# include "basic.h"
# include "scene.h"
# include<cmath>
# include<iostream>

using namespace std;

//get the specific point where the ray hit the object
Vector3 getPoint(Ray ray, double t)
{
    return(ray.start_point+ray.dir*t);
}

Plane getTriPlane(Triangle_Face tf, Vector3* vertex)
{
    Plane res;
    Vector3 normal = getFaceNormal(tf, vertex);
    res.A = normal.x;
    res.B = normal.y;
    res.C = normal.z;
    Vector3 point = vertex[tf.p1];
    res.D = -res.A*point.x - res.B*point.y - res.C*point.z;

    return res;
}

double rayPlaneIntersect(Ray ray, Plane plane)
{
    double t = 9999;
    double discriminant;
    discriminant = plane.A*ray.dir.x + plane.B*ray.dir.y + plane.C*ray.dir.z;
    if(discriminant!=0)
    {
        double numerator = plane.A*ray.start_point.x + plane.B*ray.start_point.y + plane.C*ray.start_point.z + plane.D;
        numerator *= -1.0;
        t = numerator/discriminant;
    }

    return t;
}

Vector3 barycentric(Triangle_Face tf, Vector3* vertex, Vector3 p)
{
    Vector3 res;
    Vector3 temp1 = vertex[tf.p2] - vertex[tf.p1];
    Vector3 temp2 = vertex[tf.p3] - vertex[tf.p1];
    Vector3 temp_p = p - vertex[tf.p1];
    double d11 = temp1.dotProd(temp1);
    double d12 = temp1.dotProd(temp2);
    double d22 = temp2.dotProd(temp2);
    double dp1 = temp_p.dotProd(temp1);
    double dp2 = temp_p.dotProd(temp2);
    double D = d11*d22 - d12*d12;
    double D_b = d22*dp1 - d12*dp2;
    double D_g = d11*dp2 - d12*dp1;
    res.y = D_b/D;
    res.z = D_g/D;
    res.x = 1.0-(res.y+res.z);
    return res;
}

bool isInTriangle(Vector3 p, Triangle_Face tf, Vector3* vertex)
{
    bool res = 0;
    Vector3 b_coordinate = barycentric(tf, vertex, p);
    if(b_coordinate.x>=0&&b_coordinate.y>=0&&b_coordinate.z>=0&&b_coordinate.x+b_coordinate.y+b_coordinate.z<=1.001)
        res = 1;
    return res;
}

//Decide the closest&positive ray/sphere intersection point
double closestIntersect(Sphere sphere_a, Ray ray)
{
    double B, C;
    double discriminant;
    B = 2*(ray.dir.x*(ray.start_point.x-sphere_a.center.x) + ray.dir.y*(ray.start_point.y-sphere_a.center.y) + ray.dir.z*(ray.start_point.z-sphere_a.center.z));
    C = pow(ray.start_point.x-sphere_a.center.x, 2.0)+pow(ray.start_point.y-sphere_a.center.y,2.0)+pow(ray.start_point.z-sphere_a.center.z,2.0)-pow(sphere_a.radius, 2.0);
    discriminant = B*B-4*C;
    if(discriminant < 0)
    {
        return 9999.0; //no intersection
    }
    double t1, t2, t;
    t1 = 0.5*(-B+sqrt(discriminant));
    t2 = 0.5*(-B-sqrt(discriminant));
    if(t1>0&&t2>0)
    {
        t = min(t1,t2);
    }
    else if(t1>0&&t2<=0)
    {
        t = t1;
    }
    else if(t1<=0&&t2>0)
    {
        t = t2;
    }
    else{
        t = 9999.0;
    }
    return t;
}

Color biLinInterplt(Scene_objs scene, double u, double v)
{
    Color res;
    
    Vector2 ul_indx = {floor(u*(scene.text_size.u-1)), floor(v*(scene.text_size.v-1))};
    Vector2 ur_indx = {ul_indx.u+1, ul_indx.v};
    Vector2 ll_indx = {ul_indx.u, ul_indx.v+1};
    Vector2 lr_indx = {ul_indx.u+1, ul_indx.v+1};
    
    v = v*(scene.text_size.v-1);
    u = u*(scene.text_size.u-1);
    int width = scene.text_size.u;
    int height = scene.text_size.v; 
    
    //compute intermediate point in the first row
    double rm1 = scene.texture[ul_indx.u][ul_indx.v].r *(ur_indx.u-u) + scene.texture[ur_indx.u][ur_indx.v].r *(u-ul_indx.u);   
    double gm1 = scene.texture[ul_indx.u][ul_indx.v].g *(ur_indx.u-u) + scene.texture[ur_indx.u][ur_indx.v].g *(u-ul_indx.u);
    double bm1 = scene.texture[ul_indx.u][ul_indx.v].b *(ur_indx.u-u) + scene.texture[ur_indx.u][ur_indx.v].b *(u-ul_indx.u);
    //compute intermediate point in the last row
    double rm2 = scene.texture[ll_indx.u][ll_indx.v].r *(lr_indx.u-u) + scene.texture[lr_indx.u][lr_indx.v].r *(u-ll_indx.u);
    double gm2 = scene.texture[ll_indx.u][ll_indx.v].g *(lr_indx.u-u) + scene.texture[lr_indx.u][lr_indx.v].g *(u-ll_indx.u);
    double bm2 = scene.texture[ll_indx.u][ll_indx.v].b *(lr_indx.u-u) + scene.texture[lr_indx.u][lr_indx.v].b *(u-ll_indx.u);
    //perfrom linear interpolation in column direction between intermediate points and get the final point's color
    res.r = (ll_indx.v-v)*rm1+(v-ul_indx.v)*rm2;
    res.g = (ll_indx.v-v)*gm1+(v-ul_indx.v)*gm2;
    res.b = (ll_indx.v-v)*bm1+(v-ul_indx.v)*bm2;
    return res;
}

Color nesarestNeighbor(Scene_objs scene, double u, double v)
{
    int upos = (scene.text_size.u-1)*u;
    int vpos = (scene.text_size.v-1)*v;
    Color res = scene.texture[upos][vpos];
    return res;
}

Color textMap_sph(Scene_objs scene, Sphere sphere, Vector3 intersect)
{
    double theta, phi; 
    double u, v;
    Vector3 normal = normalize(intersect-sphere.center);
    phi = acos(normal.z);   
    theta = atan(normal.y/normal.x);
    v = phi/M_PI;
    u = 0.5*(1.0+theta/M_PI);
    Color res = biLinInterplt(scene, u, v);
    return res;
}

Color textMap_Tri(Scene_objs scene, Triangle_Face tf, Vector3 intersect)
{
    Vector3 baryCtr = barycentric(tf, scene.vertex, intersect);
    double u = baryCtr.x*scene.v_coord[tf.vt1].u+baryCtr.y*scene.v_coord[tf.vt2].u+baryCtr.z*scene.v_coord[tf.vt3].u;
    double v = baryCtr.x*scene.v_coord[tf.vt1].v+baryCtr.y*scene.v_coord[tf.vt2].v+baryCtr.z*scene.v_coord[tf.vt3].v;
    Color res = biLinInterplt(scene, u, v);
    //Color res = nesarestNeighbor(scene, u, v);    //use nearest-neighbor to decide intrinsic color
    return res;
}

Ray getIncident(Ray init, Vector3 intersect)
{
    Ray I;
    Vector3 dir;
    dir = dir - init.dir;
    I.dir = dir;
    I.start_point = intersect;
    I.isInside = init.isInside;
    return I;
}

Ray getReflect(Ray inc, Vector3 normal)
{
    double angleCos = getAngleCos(inc.dir, normal);
    Vector3 A = normal*angleCos;
    Ray R;
    R.dir = A*2 - inc.dir;
    R.start_point = inc.start_point;
    R.isInside = inc.isInside;
    return R;
}

double getF0(MtlColor a)
{
    double f0;
    f0 = pow(((a.eta-1.0)/(a.eta+1.0)),2);
    return f0;
}

double getFresnel(double f0, double cos_in)
{
    double f;
    f = f0+(1.0-f0)*pow((1.0-cos_in),5);
    return f;
}

Ray getTransmit(Ray init, Vector3 normal, double indx_i, double indx_t)
{
    Ray transmit;
    transmit.start_point = init.start_point;
    Vector3 tNorm;
    tNorm = tNorm-normal;
    double angleCos = getAngleCos(init.dir, normal);
    Vector3 dir = tNorm*pow((1.0-(pow((indx_i/indx_t),2)*(1.0-angleCos*angleCos))),0.5) + (normal*angleCos-init.dir)*(indx_i/indx_t);
    transmit.dir = dir;
    transmit.isInside = 1-init.isInside;
    return transmit;
}

//Decide the color of pixels wrt illumination for spheres
Color ShadeRay(Light* lightSrc, Scene_objs scene, Camera camera, Vector3 view_init, Sphere sphere_a, Vector3 intersect)
{
    Color resColor;
    //normal vector of a intersect point can be computed by subtract sphere center from the point vector
    Vector3 N = normalize(intersect-sphere_a.center);
    Vector3 L_;
    Vector3 toV, H;
    Ray shadowRay;
    double S = 1.0;//shadow flag, update within each loop
    
    if(sphere_a.useText)
    {
        sphere_a.surf_material.r = textMap_sph(scene, sphere_a, intersect).r;
        sphere_a.surf_material.g = textMap_sph(scene, sphere_a, intersect).g;
        sphere_a.surf_material.b = textMap_sph(scene, sphere_a, intersect).b;
    }
    //ambient color
    resColor.r = sphere_a.surf_material.r*sphere_a.surf_material.ka;
    resColor.g = sphere_a.surf_material.g*sphere_a.surf_material.ka;
    resColor.b = sphere_a.surf_material.b*sphere_a.surf_material.ka;
    //Use shadow flags to decide diffuse/specular color
    //hard shadow
    int l_count = 0;
    while(lightSrc[l_count].w < 2.0)
    {
        shadowRay.start_point = intersect;
        //L for point light/directional light
        if(lightSrc[l_count].w == 1.0)
        {
            L_ = normalize(lightSrc[l_count].pos-intersect);//point light
            shadowRay.dir = L_;
            toV = normalize(view_init - intersect);
        }
        else if(lightSrc[l_count].w == 0.0)
        {
            L_ = normalize(L_ - lightSrc[l_count].pos);//directional light
            Vector3 temp;
            toV = normalize(temp - camera.viewdir);
        }
        //decide vector H
        H = normalize(L_ + toV);
        //decide S for the current shadowRay
        int sph_count = 0;
        while(scene.sphere[sph_count].radius > 0)
        {
            if(sphere_a == scene.sphere[sph_count])
            {
                sph_count++;
                continue;
            }
            else
            {
                double t_ = closestIntersect(scene.sphere[sph_count], shadowRay);//check if the intersect is visible to the current light source
                if(t_ > 0&&t_ < 9999.0) //exist an intersection along the way
                {
                    S = 0.0;
                    break;
                }
            }
            sph_count++;
        }
        int s_count = 1;
        while(scene.triangle[s_count].p1 > 0)
        {
            Plane currPlane = getTriPlane(scene.triangle[s_count], scene.vertex);
            double t_ = rayPlaneIntersect(shadowRay, currPlane);//check if the intersect is visible to the current light source
            if(t_ > 0&&t_ < 9999.0) //exist an intersection along the way
            {
                Vector3 temp = getPoint(shadowRay, t_);
                if(isInTriangle(temp, scene.triangle[s_count], scene.vertex)==1);
                {
                    S = 0.0;
                    break;
                }                   
            }
            s_count++;
        }
        if(lightSrc[l_count].isAtt)
        {
            double d = (lightSrc[l_count].pos-intersect).lenVec();
            double F = getAttF(lightSrc[l_count], d);
            S *= F;
        }
        //Compute color based on Phong illumination model
        resColor.r += S*(max((N.dotProd(L_)),0.0)*sphere_a.surf_material.r*sphere_a.surf_material.kd + pow(max((N.dotProd(H)),0.0),sphere_a.surf_material.exp_n)*sphere_a.surf_material.ks*sphere_a.surf_material.Osr);
        resColor.g += S*(max((N.dotProd(L_)),0.0)*sphere_a.surf_material.g*sphere_a.surf_material.kd + pow(max((N.dotProd(H)),0.0),sphere_a.surf_material.exp_n)*sphere_a.surf_material.ks*sphere_a.surf_material.Osg);
        resColor.b += S*(max((N.dotProd(L_)),0.0)*sphere_a.surf_material.b*sphere_a.surf_material.kd + pow(max((N.dotProd(H)),0.0),sphere_a.surf_material.exp_n)*sphere_a.surf_material.ks*sphere_a.surf_material.Osb);
        l_count++;
        S = 1.0;
    }
    
    //clamp the color value so that it doesn't exceed 1.0
    /*
    resColor.r = min(resColor.r, 1.0);
    resColor.g = min(resColor.g, 1.0);
    resColor.b = min(resColor.b, 1.0);*/
    return resColor;
}

Color ShadeRayTri(Light* lightSrc, Scene_objs scene, Camera camera, Triangle_Face tf, Vector3 intersect)
{
    Color resColor;
    //normal vector of a intersect point can be computed by subtract sphere center from the point vector
    Vector3 N = getFaceNormal(tf, scene.vertex);
    //Vector3 N;
    Vector3 L_;
    Vector3 toV, H;
    Vector3 baryCtr = barycentric(tf, scene.vertex, intersect);
    Ray shadowRay;
    double S = 1.0;//shadow flag, update within each loop
    if(tf.useText)
    {
        tf.material.r = textMap_Tri(scene, tf, intersect).r;
        tf.material.g = textMap_Tri(scene, tf, intersect).g;
        tf.material.b = textMap_Tri(scene, tf, intersect).b;
    }
    //ambient color
    resColor.r = tf.material.r*tf.material.ka;
    resColor.g = tf.material.g*tf.material.ka;
    resColor.b = tf.material.b*tf.material.ka;
    //Use shadow flags to decide diffuse/specular color
    //hard shadow
    int l_count = 0;
    while(lightSrc[l_count].w < 2.0)
    {
        shadowRay.start_point = intersect;
        if(tf.isFlat == false)
        {
            int v1 = tf.normal1;
            int v2 = tf.normal2;
            int v3 = tf.normal3;
            N = getFaceNormal(tf, scene.vertex);
            Vector3 zero;
            if(scene.v_normal[v1] != zero)
            {
                Vector3 w1 = normalize(scene.v_normal[v1]);
                Vector3 w2 = normalize(scene.v_normal[v2]);
                Vector3 w3 = normalize(scene.v_normal[v3]);
                N = normalize(w1*baryCtr.x+w2*baryCtr.y+w3*baryCtr.z);
            }
        }    
        //L for point light/directional light
        if(lightSrc[l_count].w == 1.0)
        {                    
            L_ = normalize(lightSrc[l_count].pos-intersect);//point light
            shadowRay.dir = L_;
            toV = normalize(camera.eye - intersect);
        }
        else if(lightSrc[l_count].w == 0.0)
        {
            L_ = normalize(L_ - lightSrc[l_count].pos);//directional light
            Vector3 temp;
            toV = normalize(temp - camera.viewdir);
        }
        //decide vector H
        H = normalize(L_ + toV);
        //decide S for the current shadowRay
        /*
        int sph_count = 0;       
        while(scene.sphere[sph_count].radius>0)
        {
            double t2 = closestIntersect(scene.sphere[sph_count], shadowRay);
            if(t2 > 0&&t2 < 9999.0) //exist an intersection along the way
            {
                S = 0.0;
                break;                 
            }
            sph_count++;
        }
        int s_count = 1;
        while(scene.triangle[s_count].p1 > 0)
        {
            if(tf == scene.triangle[s_count])
            {
                s_count++;
                continue;
            }
            else
            {
                Plane currPlane = getTriPlane(scene.triangle[s_count], scene.vertex);
                double t_ = rayPlaneIntersect(shadowRay, currPlane);//check if the intersect is visible to the current light source
                if(t_ > 0&&t_ < 9999.0) //exist an intersection along the way
                {
                    Vector3 temp = getPoint(shadowRay, t_);
                    if(isInTriangle(temp, scene.triangle[s_count], scene.vertex)==1);
                    {
                        S = 0.0;
                        break;
                    }                   
                }
            }
            s_count++;
        }*/
        if(lightSrc[l_count].isAtt)
        {
            double d = (lightSrc[l_count].pos-intersect).lenVec();
            double F = getAttF(lightSrc[l_count], d);
            S *= F;
        }
        //Compute color based on Phong illumination model
        resColor.r += S*(max((N.dotProd(L_)),0.0)*tf.material.r*tf.material.kd + pow(max((N.dotProd(H)),0.0),tf.material.exp_n)*tf.material.ks*tf.material.Osr);
        resColor.g += S*(max((N.dotProd(L_)),0.0)*tf.material.g*tf.material.kd + pow(max((N.dotProd(H)),0.0),tf.material.exp_n)*tf.material.ks*tf.material.Osg);
        resColor.b += S*(max((N.dotProd(L_)),0.0)*tf.material.b*tf.material.kd + pow(max((N.dotProd(H)),0.0),tf.material.exp_n)*tf.material.ks*tf.material.Osb);
        l_count++;
        S = 1.0;
    }
    //clamp the color value so that it doesn't exceed 1.0
    /*
    resColor.r = min(resColor.r, 1.0);
    resColor.g = min(resColor.g, 1.0);
    resColor.b = min(resColor.b, 1.0);*/
    return resColor;
}

Color TraceSingle(Light* lightSrc, Scene_objs scene, Camera camera, int maxDepth, int depth, Ray cur_ray)
{
    if(depth > maxDepth)
    {
        return scene.background;
    }
    Ray I, R, T;
    Vector3 curr_inter,new_inter, trans_inter, surfN, view_init;
    Color resColor;
    double F;
    double f0, cosin;
    int hit_type = 0;//1 for hit sphere, 2 for hit triangle

    int k=0; //For the ray passes one specific pixel, check all spheres in the scene for intersection
    double intersect = 9999.0; //Scale parameter t in the ray function at intersection point. Updates through loop to find the minimal
    int hit_k = 0;
    while(scene.sphere[k].radius > 0)
    {
        double thisIntersect = closestIntersect(scene.sphere[k], cur_ray);
        if(thisIntersect < intersect)
        {
            intersect = thisIntersect;
            hit_k = k;
            hit_type = 1;            
        }
        k++;
    }
    
    int num_tri = 1;
    double intersect2 = 9999.0;
    int hit_t = 0;
    while(scene.triangle[num_tri].p1>0)
    {
        Plane triPlane = getTriPlane(scene.triangle[num_tri], scene.vertex);
        double thisIntersect = rayPlaneIntersect(cur_ray, triPlane);
        if(thisIntersect < intersect2)
        {
            Vector3 intPoint = getPoint(cur_ray, thisIntersect);
            if(isInTriangle(intPoint, scene.triangle[num_tri], scene.vertex)==1)
            {
                intersect2 = thisIntersect;
                hit_t = num_tri;
            }           
        }
        num_tri++;
    }    
    if(intersect2 < intersect)
    {
        hit_type = 2;
    }
    
    if(hit_type == 0)
    {
        resColor = scene.background;
        return resColor;
    }

    if(hit_type == 1)
    {
        new_inter = getPoint(cur_ray, intersect);
        I = getIncident(cur_ray, new_inter);
        surfN = normalize(new_inter-scene.sphere[hit_k].center); //surface normal of intersection(only currnt sphere)
        if(I.dir.dotProd(surfN) < 0)
        {
            Vector3 temp;
            surfN = temp - surfN;
        }
        cosin = getAngleCos(I.dir, surfN);       
        f0 = getF0(scene.sphere[hit_k].surf_material);
        F = getFresnel(f0, cosin);
        
        Color preColor = ShadeRay(lightSrc, scene, camera, view_init, scene.sphere[hit_k], new_inter);
        resColor.r = preColor.r;
        resColor.g = preColor.g;
        resColor.b = preColor.b;

        R = getReflect(I, surfN);
        R.start_point = R.start_point+surfN*0.001;
        Color refIllum = TraceSingle(lightSrc, scene, camera, maxDepth, depth+1, R);
        resColor.r += refIllum.r*F;
        resColor.g += refIllum.g*F;
        resColor.b += refIllum.b*F;
               
        double trans_factor = (1.0-F)*(1.0-scene.sphere[hit_k].surf_material.alpha);
        if(trans_factor != 0)
        {
            double inci_indx = scene.sphere[hit_k].surf_material.eta;
            double tran_indx = 1.0;            
            if(I.isInside == 0)
            {
                inci_indx = 1.0;
                tran_indx = scene.sphere[hit_k].surf_material.eta;
            }
            T = getTransmit(I, surfN, inci_indx, tran_indx);
            T.start_point = T.start_point-surfN*0.001;
            Color transIllum = TraceSingle(lightSrc, scene, camera, maxDepth, depth+1, T);
            resColor.r += transIllum.r*trans_factor;
            resColor.g += transIllum.g*trans_factor;
            resColor.b += transIllum.b*trans_factor;
        }       
        
    }
    if(hit_type == 2)
    {
        new_inter = getPoint(cur_ray, intersect2);
        I = getIncident(cur_ray, new_inter);
        surfN = getFaceNormal(scene.triangle[hit_t], scene.vertex);
        if(I.dir.dotProd(surfN) < 0)
        {
            Vector3 temp;
            surfN = temp - surfN;
        }
        R = getReflect(I, surfN);
        f0 = getF0(scene.triangle[hit_t].material);
        cosin = getAngleCos(I.dir, surfN);
        F = getFresnel(f0, cosin);
        
        Color preColor = ShadeRayTri(lightSrc, scene, camera, scene.triangle[hit_t], new_inter);
        resColor.r = preColor.r;
        resColor.g = preColor.g;
        resColor.b = preColor.b;

        new_inter = R.start_point+surfN*0.001;
        R.start_point = new_inter;
        Color refIllum = TraceSingle(lightSrc, scene, camera, maxDepth, depth+1, R);
        resColor.r += refIllum.r*F;
        resColor.g += refIllum.g*F;
        resColor.b += refIllum.b*F;

        double trans_factor = (1.0-F)*(1.0-scene.triangle[hit_t].material.alpha);
        if(trans_factor != 0)
        {
            double inci_indx = 1.0;
            double tran_indx = scene.triangle[hit_t].material.eta;
            if(I.isInside == 1)
            {
                inci_indx = scene.triangle[hit_t].material.eta;
                tran_indx = 1.0;
            }
            T = getTransmit(I, surfN, inci_indx, tran_indx);
            T.start_point = T.start_point-surfN*0.001;
            Color transIllum = TraceSingle(lightSrc, scene, camera, maxDepth, depth+1, T);
            resColor.r += transIllum.r*trans_factor;
            resColor.g += transIllum.g*trans_factor;
            resColor.b += transIllum.b*trans_factor;
        }       
    }
    return resColor;
}

//Decide color for each pixel in the image
//For each pixel in the image, generate a ray that starts from the eye and goes across the pixel.
void RayTracer(Color** renderPic, int width, int height, Color bkg_color, Camera camera, Light* lightSrc, Scene_objs scene, int maxLevel, Vector3 ul, Vector3 dh, Vector3 dv)
{
    int maxDepth = maxLevel;
    Ray curr_ray;
    curr_ray.start_point = camera.eye;
    for(int j=0; j < height; j++)
    {
        for(int i=0; i < width; i++)
        {
            renderPic[i][j] = bkg_color;
            
            Vector3 pos_v = ul+dh*(i+0.5)+dv*(j+0.5);
            curr_ray.dir = normalize(pos_v-curr_ray.start_point);

            Color illumiColor = TraceSingle(lightSrc, scene, camera, maxDepth, 0, curr_ray);
            renderPic[i][j].r = illumiColor.r;
            renderPic[i][j].g = illumiColor.g;
            renderPic[i][j].b = illumiColor.b;
        }
    }
}

void clamp(Color** renderPic, int width, int height)
{
    double maxR, maxG, maxB;
    maxR = renderPic[0][0].r;
    maxG = renderPic[0][0].g;
    maxB = renderPic[0][0].b;
    
    for(int j=0; j < height; j++)
    {
        for(int i=0; i < width; i++)
        {
            if(maxR < renderPic[i][j].r)
            {
                maxR = renderPic[i][j].r;
            }
            if(maxG < renderPic[i][j].g)
            {
                maxG = renderPic[i][j].g;
            }
            if(maxB < renderPic[i][j].b)
            {
                maxB = renderPic[i][j].b;
            }            
        }
    }
    
    for(int j=0; j < height; j++)
    {
        for(int i=0; i < width; i++)
        {
            renderPic[i][j].r /= maxR;
            renderPic[i][j].g /= maxG;
            renderPic[i][j].b /= maxB;
            
            renderPic[i][j].r = floor(renderPic[i][j].r*255.0);
            renderPic[i][j].g = floor(renderPic[i][j].g*255.0);
            renderPic[i][j].b = floor(renderPic[i][j].b*255.0);           
        }
    }
}