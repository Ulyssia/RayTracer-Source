# define _USE_MATH_DEFINES
# include <iostream>
# include <fstream>
# include <string>
# include <cmath>
# include "funcs.h"
# include "basic.h"
# include "scene.h"

using namespace std;

int main(int args, char* argv[])
{
    std::string type = "P3";
    std::string inputFile, key;
    Color** myRender; //the image pixel array

    Vector3 u, v, n;//define viewing window
    double asp_ratio;
    double view_dist = 10.0; //viewing distance
    double vw, vh; //width and height of viewing window
    int width = 0, height = 0;

    Color bkg_color; 
    MtlColor mtl_color[10]; //all params for illumination equation
    Light lightSource[10];
    Camera myCamera;
    Scene_objs myScene;
    Sphere mySphere[10]; //allow rendering maximum of 10 spheres in the scene
    Vector3 vertex_list[200];
    Vector3 vnormal_list[200];
    Vector2 vt_list[200];
    Triangle_Face triangle[100];

    int sphere_count = 0;
    int vertex_num = 1;
    int vn_num = 1;
    int vt_num = 1;
    int tri_num = 1;
    int material_count = 0;
    int light_count = 0;
    int maxLevel = 0;
    //Open the description file
    inputFile = argv[1];
    maxLevel = stoi(argv[2]);

    int texture_w=0, texture_h=0;
    Color** textureMap;
    //bool useTexture = false;

    //Read in all parameters from input file
    ifstream desFile(inputFile.append(".txt"));
    if(!desFile)
    {
        cout << "Fail opening input description file!" << endl;
    }
    while(desFile >> key)
    {
        if(key == "imsize")
        {
            desFile >> width >> height;
            
            if(width == 0 || height == 0)
            {
                cout << "Incorrect or missing information of image size! Check input file." << endl;
                break;
            }
            //define a pixel array of height*width
            myRender = new Color*[width];
            for(int w = 0; w<width; w++)
            {
                myRender[w] = new Color[height];
            }
        }
        if(key == "eye")
        {
            desFile >> myCamera.eye.x >> myCamera.eye.y >> myCamera.eye.z;
        }
        if(key == "viewdir")
        {
            desFile >> myCamera.viewdir.x >> myCamera.viewdir.y >> myCamera.viewdir.z;
        }
        if(key == "updir")
        {
            desFile >> myCamera.updir.x >> myCamera.updir.y >> myCamera.updir.z;
        }
        if(key == "vfov")
        {
            desFile >> myCamera.vfov;
        }
        if(key == "light" || key == "attlight")
        {
            desFile >> lightSource[light_count].pos.x >> lightSource[light_count].pos.y >> lightSource[light_count].pos.z >> lightSource[light_count].w >> lightSource[light_count].light_color.r >> lightSource[light_count].light_color.g >> lightSource[light_count].light_color.b;
            if(key == "attlight")
            {
                lightSource[light_count].isAtt = 1;
                desFile >> lightSource[light_count].attParam.x >> lightSource[light_count].attParam.y >> lightSource[light_count].attParam.z;
            }
            light_count++;
        }
        if(key == "bkgcolor")
        {
            desFile >> bkg_color.r >> bkg_color.g >> bkg_color.b;
        }
        if(key == "mtlcolor")
        {
            desFile >> mtl_color[material_count].r >> mtl_color[material_count].g >> mtl_color[material_count].b >> mtl_color[material_count].Osr >> mtl_color[material_count].Osg >> mtl_color[material_count].Osb >> mtl_color[material_count].ka >> mtl_color[material_count].kd >> mtl_color[material_count].ks >> mtl_color[material_count].exp_n >> mtl_color[material_count].alpha >> mtl_color[material_count].eta;
            material_count++;
        }
        if(key == "sphere")
        {
            desFile >> mySphere[sphere_count].center.x >> mySphere[sphere_count].center.y >> mySphere[sphere_count].center.z;
            desFile >> mySphere[sphere_count].radius;
            mySphere[sphere_count].surf_material = mtl_color[material_count-1];
            mySphere[sphere_count].surf_material.r = mtl_color[material_count-1].r;
            mySphere[sphere_count].surf_material.g = mtl_color[material_count-1].g;
            mySphere[sphere_count].surf_material.b = mtl_color[material_count-1].b;
            if(texture_w>0&&texture_h>0)
            {
                mySphere[sphere_count].useText = true;
            }
            sphere_count ++;
        }
        if(key == "v")
        {
            desFile >> vertex_list[vertex_num].x >> vertex_list[vertex_num].y >> vertex_list[vertex_num].z;
            vertex_num ++;
        }
        if(key == "vn")
        {
            desFile >> vnormal_list[vn_num].x >> vnormal_list[vn_num].y >> vnormal_list[vn_num].z;
            vn_num++;
        }
        if(key == "f")
        {
            string v_c[3];
            for(int i = 0; i < 3; i++)
            {
                desFile >> v_c[i];
            }
            if(texture_w>0&&texture_h>0)
            {
                triangle[tri_num].useText = true;
            }
            triangle[tri_num].p1 = stoi(string(1,v_c[0][0]));
            triangle[tri_num].p2 = stoi(string(1,v_c[1][0]));          
            triangle[tri_num].p3 = stoi(string(1,v_c[2][0]));
            if(v_c[0].length()>=3)
            {
                if(triangle[tri_num].useText)
                {
                    triangle[tri_num].vt1 = stoi(string(1,v_c[0][2]));
                    triangle[tri_num].vt2 = stoi(string(1,v_c[1][2]));
                    triangle[tri_num].vt3 = stoi(string(1,v_c[2][2]));
                    if(v_c[0].length()==5)
                    {
                        triangle[tri_num].normal1 = stoi(string(1,v_c[0][4]));
                        triangle[tri_num].normal2 = stoi(string(1,v_c[1][4]));
                        triangle[tri_num].normal3 = stoi(string(1,v_c[2][4]));
                        triangle[tri_num].isFlat = false;
                    }
                }
                else
                {
                    triangle[tri_num].normal1 = stoi(string(1,v_c[0][2]));
                    triangle[tri_num].normal2 = stoi(string(1,v_c[1][2]));
                    triangle[tri_num].normal3 = stoi(string(1,v_c[2][2]));             
                    triangle[tri_num].isFlat = false;
                }
            }           
            //desFile >> triangle[tri_num].p1 >> triangle[tri_num].p2 >> triangle[tri_num].p3;
            triangle[tri_num].material = mtl_color[material_count-1];
            triangle[tri_num].material.r = mtl_color[material_count-1].r;
            triangle[tri_num].material.g = mtl_color[material_count-1].g;
            triangle[tri_num].material.b = mtl_color[material_count-1].b;
            tri_num ++;
        }
        if(key == "texture")
        {
            string textureName;
            desFile >> textureName;
            ifstream textureFile(textureName);
            int line = 1;
            string temp1, temp2;
            textureFile >> temp1 >> texture_w >> texture_h >> temp2;
            textureMap = new Color*[texture_w];
            for(int w = 0; w<texture_w; w++)
            {
                textureMap[w] = new Color[texture_h];
            }
            for(int i = 0; i < texture_h; i++)
            {
                for(int j = 0; j < texture_w; j++)
                {
                    textureFile >> textureMap[j][i].r >> textureMap[j][i].g >> textureMap[j][i].b;
                    textureMap[j][i].r /= 255.0;
                    textureMap[j][i].g /= 255.0;
                    textureMap[j][i].b /= 255.0;
                }
            }
            textureFile.close();
        }
        if(key == "vt")
        {
            desFile >> vt_list[vt_num].u >> vt_list[vt_num].v;
            vt_num++;
        }
        if(key == "viewdist")
        {
            desFile >> view_dist;
        }
    }                                                                                                   
    desFile.close();

    //Defining the viewing window  
    u = myCamera.viewdir.crossProd(myCamera.updir);
    u = normalize(u);
    v = u.crossProd(myCamera.viewdir);
    v = normalize(v);
    n = normalize(myCamera.viewdir);
    asp_ratio = (double)1.00*width/height;
    //view_dist = 12;
    vh = 2.0*view_dist*tan(myCamera.vfov*M_PI/360.0);
    vw = 1.0*vh*asp_ratio;

    Vector3 ul, ur, ll, lr; //viewing window corners
    ul = myCamera.eye + n*view_dist - u*(vw/2.0) + v*(vh/2.0);
    ur = myCamera.eye + n*view_dist + u*(vw/2.0) + v*(vh/2.0);
    ll = myCamera.eye + n*view_dist - u*(vw/2.0) - v*(vh/2.0);
    lr = myCamera.eye + n*view_dist + u*(vw/2.0) - v*(vh/2.0);

    Vector3 delta_v, delta_h; //Steps used when go across the image pixels
    delta_h = (ur-ul)*(1.0/width);
    delta_v = (ll-ul)*(1.0/height);

    myScene.sphere = mySphere;
    myScene.triangle = triangle;
    myScene.vertex = vertex_list;
    myScene.v_normal = vnormal_list;
    myScene.v_coord = vt_list;
    myScene.text_size = {texture_w, texture_h};
    myScene.texture = textureMap;
    myScene.background = bkg_color;

    //cout << myScene.triangle[2].p1 << myScene.triangle[2].p2 << myScene.triangle[2].p3 << endl;
    //Core step of performing raytracing
    RayTracer(myRender, width, height, bkg_color, myCamera, lightSource, myScene, maxLevel, ul, delta_h, delta_v);
    clamp(myRender,width,height);
    
    //Generate a new PPM file 
    ofstream MyPPM;
    MyPPM.open("MyImage.ppm", ios::out);
    if(!MyPPM)
    {
        cout << "Fail opening new ppm file!" << endl;
    }
    //Write the file head
    MyPPM << type << endl;
    MyPPM << "# PPM file type is P3 by default" << endl;
    MyPPM << width << " " << height << endl;
    MyPPM << 255 << endl;
    //run a double-loop to copy pixel color from array to PPM
    for(int h = 0; h < height; h++)
    {
        for(int w = 0; w < width; w++)
        {                 
            if(w > 0 && w%5 == 0)
            {
                MyPPM << "\n";
            }             
            MyPPM << myRender[w][h].r << " " << myRender[w][h].g << " " << myRender[w][h].b << " ";
        }
        MyPPM << "\n";
    }
    MyPPM.close();

    return 0;
}

