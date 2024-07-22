/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Srija Madarapu>
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif


#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>

using std::sqrt;
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0
#define PI 3.14159265

unsigned char buffer[HEIGHT][WIDTH][3];

double a = double(WIDTH) / double(HEIGHT);
double rad = fov*PI /180;


struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};


struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

struct intersect
{
    glm::vec3 position;
    glm::vec3 normal;
    float t= std::numeric_limits<float>::max();
    int sphereID;
    int triID;
    glm::vec3 triabg;
    int shape; //0 sphere 1 triangle
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

class ray {
public:
    ray() {}
    ray(glm::vec3& origin, glm::vec3& direction) : orig(origin), dir(direction) {}

    glm::vec3 origin() { return orig; }
    glm::vec3 direction() { return dir; }

    glm::vec3 intersect_at(float t) {
        return orig + t * dir;
    }
    glm::vec3 orig;
    glm::vec3 dir;
};

bool hit_sphere(glm::vec3 center, double radius, ray& r,intersect& intersect_p) {
    glm::vec3 oc = center-r.origin();
    double a = glm::dot(r.direction(), r.direction()) ;
    double h = glm::dot(oc, r.direction());
    double c = glm::dot(oc, oc) - radius * radius;
    double discriminant = h * h -  a * c;
    float t = 0;
    if (discriminant < 0) {
        return false;
    }
    else {
        t = glm::min((h - sqrt(discriminant)) / (a), (h + sqrt(discriminant)) / (a));
        glm::vec3 inter = r.intersect_at(t);
        intersect_p.position = inter;
        glm::vec3 n = glm::normalize((inter - center));
        intersect_p.normal = n;
        intersect_p.t = t;
        return true;
    }  
}
float area(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2)
{
    glm::vec3 e1 = v1 - v0;
    glm::vec3 e2 = v2 - v0;
    glm::vec3 crossProduct = glm::cross(e1, e2);
    return (glm::sqrt(glm::dot(crossProduct, crossProduct)));
}
glm::vec3 tri_point(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2,glm::vec3 p)
{
    float ABC = area(v0, v1, v2);

    float PBC = area(p, v1, v2);
    float APC = area(v0, p, v2);
    float ABP = area(v0, v1, p);

    return (glm::vec3(PBC / ABC, APC / ABC, ABP / ABC));
}

double hit_triangle(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 n0, glm::vec3 n1, glm::vec3 n2, ray& r,intersect& I) {
    glm::vec3 e1 = v1 - v0;
    glm::vec3 e2 = v2 - v0;
    glm::vec3 pvec = glm::cross(r.dir, e2);
    double det = glm::dot(e1, pvec);

    // Check if the ray and triangle are parallel
    if (glm::abs(det) < 0.0001) {
        return false;
    }

    double invDet = 1.0 / det;
    glm::vec3 tvec = r.orig - v0;
    double u = glm::dot(tvec, pvec) * invDet;
    if (u < 0.0 || u > 1.0) {
        return false;
    }

    glm::vec3 qvec = glm::cross(tvec, e1);
    double v = glm::dot(r.dir, qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }

    double t = glm::dot(e2, qvec) * invDet;
    if (t > 0.0000000000001) {
        glm::vec3 inter = r.intersect_at(t);
        I.position = inter;
        glm::vec3 n = glm::normalize(glm::cross(e1,e2));
        I.normal = n;
        I.t = t;
        return true;
    }
    else
    {
        return false;
    }
}
glm::vec3 barycentric_coordinates(intersect I,int n)
{
    Triangle tri = triangles[I.triID];
    glm::vec3 a;
    switch (n)
    {
    case 0:
        a=glm::vec3(tri.v[0].normal[0] * I.triabg[0] + tri.v[1].normal[0] * I.triabg[1] + tri.v[2].normal[0] * I.triabg[2],
        tri.v[0].normal[1] * I.triabg[0] + tri.v[1].normal[1] * I.triabg[1] + tri.v[2].normal[1] * I.triabg[2],
        tri.v[0].normal[2] * I.triabg[0] + tri.v[1].normal[2] * I.triabg[1] + tri.v[2].normal[2] * I.triabg[2]);
        break;
    case 1:
        a = glm::vec3(tri.v[0].color_diffuse[0] * I.triabg[0] + tri.v[1].color_diffuse[0] * I.triabg[1] + tri.v[2].color_diffuse[0] * I.triabg[2],
            tri.v[0].color_diffuse[1] * I.triabg[0] + tri.v[1].color_diffuse[1] * I.triabg[1] + tri.v[2].color_diffuse[1] * I.triabg[2],
            tri.v[0].color_diffuse[2] * I.triabg[0] + tri.v[1].color_diffuse[2] * I.triabg[1] + tri.v[2].color_diffuse[2] * I.triabg[2]);
        break;
    case 2:
        a = glm::vec3(tri.v[0].color_specular[0] * I.triabg[0] + tri.v[1].color_specular[0] * I.triabg[1] + tri.v[2].color_specular[0] * I.triabg[2],
            tri.v[0].color_specular[1] * I.triabg[0] + tri.v[1].color_specular[1] * I.triabg[1] + tri.v[2].color_specular[1] * I.triabg[2],
            tri.v[0].color_specular[2] * I.triabg[0] + tri.v[1].color_specular[2] * I.triabg[1] + tri.v[2].color_specular[2] * I.triabg[2]);
        break;
    default:
        break;
    }
    return a;

}
float shin(intersect I)
{
    Triangle tri = triangles[I.triID];
    return (tri.v[0].shininess* I.triabg[0] + tri.v[1].shininess * I.triabg[1] + tri.v[2].shininess * I.triabg[2]);
}
using color = glm::vec3;
void phong_shading(intersect intersection_p, Light l, glm::vec3 cam, color& col)
{
    glm::vec3 L(l.position[0], l.position[1], l.position[2]);
    L =glm::normalize(L  - intersection_p.position);

    glm::vec3 E(glm::normalize(cam - intersection_p.position));
    glm::vec3 N;
    glm::vec3 kd;
    glm::vec3 ks;
    float n;
    if (intersection_p.shape == 0)
    {
        N=glm::vec3((intersection_p.normal));
        if (glm::dot(N, L) < 0)
        {
            N = glm::vec3(0, 0, 0);
        }
        N = glm::normalize(N);
        kd=glm::vec3(spheres[intersection_p.sphereID].color_diffuse[0], spheres[intersection_p.sphereID].color_diffuse[1], spheres[intersection_p.sphereID].color_diffuse[2]);
        ks=glm::vec3(spheres[intersection_p.sphereID].color_specular[0], spheres[intersection_p.sphereID].color_specular[1], spheres[intersection_p.sphereID].color_specular[2]);
        n = spheres[intersection_p.sphereID].shininess;
    }
    else
    {
        N= barycentric_coordinates(intersection_p,0);
        kd = barycentric_coordinates(intersection_p,1);
        ks = barycentric_coordinates(intersection_p,2);
        n = shin(intersection_p);
        if (glm::dot(N, L) < 0)
        {
            N = glm::vec3(0, 0, 0);
        }
        N = glm::normalize(N);
        
    }
    float d = glm::max<float>(glm::dot(N, L),0);
    d = glm::min<float>(d, 1);

    //glm::vec3 H = L + E;
    //H = glm::normalize(H);

    glm::vec3 R = 2 * glm::dot(N, L) * N - L;
    R = glm::normalize(R);

    if (glm::dot(R, E) < 0)
    {
        R = glm::vec3(0, 0, 0);
    }

    float s = glm::max<float>(glm::dot(R, E), 0);
    s = glm::min<float>(s, 1);
    s = pow(s, n);

    glm::vec3 le(l.color[0], l.color[1], l.color[2]);
    col = col+d * kd*le+ s * ks*le;
}
color ray_color( ray& r) 
{
    intersect intersect_p;
    bool t = false;
    ////sphere
    for (int c = 0; c < num_spheres; c++)
    {
        intersect I;
        glm::vec3 p(spheres[c].position[0], spheres[c].position[1], spheres[c].position[2]);
        if (hit_sphere(p, spheres[c].radius, r, I) )
        {
            if (intersect_p.t > I.t)
            {
                intersect_p = I;
                intersect_p.sphereID = c;
                t = true;
                intersect_p.shape = 0;
            }
        }
    }
    
   //triangle
   for (int c = 0; c < num_triangles; c++)
    {
        intersect I;
        Triangle tri = triangles[c];
        Vertex vert0 = (tri.v[0]);
        Vertex vert1 = (tri.v[1]);
        Vertex vert2 = (tri.v[2]);
        glm::vec3 v0(vert0.position[0], vert0.position[1], vert0.position[2]);
        glm::vec3 v1(vert1.position[0], vert1.position[1], vert1.position[2]);
        glm::vec3 v2(vert2.position[0], vert2.position[1], vert2.position[2]);

        glm::vec3 n0(vert0.normal[0], vert0.normal[1], vert0.normal[2]);
        glm::vec3 n1(vert1.normal[0], vert1.normal[1], vert1.normal[2]);
        glm::vec3 n2(vert2.normal[0], vert2.normal[1], vert2.normal[2]);

        if (hit_triangle(v0, v1, v2, n0, n1, n2, r,I))
        {   
            if (intersect_p.t > I.t)
            {
                intersect_p = I;
                intersect_p.triID = c;
                intersect_p.triabg = tri_point(v0, v1, v2, intersect_p.position);
                intersect_p.normal = barycentric_coordinates(intersect_p, 0);
                t = true;
                intersect_p.shape = 1;
            }
        }
    }
    color col(0,0,0);
    if (!t)
    {
        col = glm::vec3(1, 1, 1);
        
    }
    else
    {
        col[0] += ambient_light[0];
        col[1] += ambient_light[1];
        col[2] += ambient_light[2];
        for (int c = 0; c < num_lights; c++)
        {
            intersect light_inter;  
            glm::vec3 L(lights[c].position[0], lights[c].position[1], lights[c].position[2]);
            float light_dis = glm::sqrt(glm::dot(L-intersect_p.position, L - intersect_p.position));
            glm::vec3 d(-glm::normalize(L-intersect_p.position));
            glm::vec3 o = L;
            
            ray reflect(o, d);
            bool t1 = false;
            for (int c = 0; c < num_spheres; c++)
            {
                intersect I;
                glm::vec3 p(spheres[c].position[0], spheres[c].position[1], spheres[c].position[2]);
                if (hit_sphere(p, spheres[c].radius, reflect, I))
                { 
                    
                    if (light_inter.t > I.t)
                    {
                        light_inter = I;
                        light_inter.sphereID = c;
                        t1 = true;
                        light_inter.shape = 0;
                    }
                }
            }
           
            for (int c = 0; c < num_triangles; c++)
            {
                intersect I;
                Triangle tri = triangles[c];
                Vertex vert0 = (tri.v[0]);
                Vertex vert1 = (tri.v[1]);
                Vertex vert2 = (tri.v[2]);
                glm::vec3 v0(vert0.position[0], vert0.position[1], vert0.position[2]);
                glm::vec3 v1(vert1.position[0], vert1.position[1], vert1.position[2]);
                glm::vec3 v2(vert2.position[0], vert2.position[1], vert2.position[2]);

                glm::vec3 n0(vert0.normal[0], vert0.normal[1], vert0.normal[2]);
                glm::vec3 n1(vert1.normal[0], vert1.normal[1], vert1.normal[2]);
                glm::vec3 n2(vert2.normal[0], vert2.normal[1], vert2.normal[2]);
                
                if (hit_triangle(v0, v1, v2, n0, n1, n2, reflect, I))
                {
                   
                    if (light_inter.t > I.t)
                    {   
                        light_inter = I;
                        light_inter.triID = c;
                        light_inter.triabg = tri_point(v0, v1, v2, light_inter.position);
                        light_inter.normal = barycentric_coordinates(light_inter, 0);
                        t1 = true;
                        light_inter.shape = 1;
                    }
                }
            }

            if (t1 && (glm::sqrt(glm::dot(light_inter.position-o,light_inter.position-o) )  <light_dis) && glm::sqrt(glm::distance2(intersect_p.position, light_inter.position))>0.001)
            {
                continue;
            }
            else
            {
                phong_shading(intersect_p, lights[c], r.orig, col);
            }

        }
        
    }
    col[0] = glm::max(col[0], 0.0f);
    col[0] = glm::min(col[0], 1.0f);

    col[1] = glm::max(col[1], 0.0f);
    col[1] = glm::min(col[1], 1.0f);

    col[2] = glm::max(col[2], 0.0f);
    col[2] = glm::min(col[2], 1.0f);
    return col;
}


//MODIFY THIS FUNCTION
void draw_scene()
{
  //a simple test output
  for(int i=0; i<WIDTH; i++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(int j=0; j<HEIGHT; j++)
    {
        double x = (2 * ((i + 0.5) / double(WIDTH)) - 1);
        double y = (2 * ((j + 0.5) / double(HEIGHT)) - 1);
        glm::vec3 origin(0, 0, 0);
        glm::vec3 dir(x*a*tan(rad/2),y*tan(rad/2),-1);
        dir = glm::normalize(dir);
        ray r(origin,dir);
        glm::vec3 col = ray_color(r);
        plot_pixel(i, j, col[0] * 255, col[1] * 255, col[2] * 255); 
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

