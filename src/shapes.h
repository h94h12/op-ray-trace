/*
 * shapes
*/


#include <iostream>
#include <math.h>
#include <vector>
#include <cfloat>
#include "primitives.h"

#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2

#define LARGE_NUM 1000000
#define NO_INTERSECTION FLT_MAX

class BRDF{
    public:
        Color kd, ks, ke;
        float kr; 
        int sp; 
    
    BRDF(); 
    BRDF(Color, Color, Color, float, int);
    BRDF clone();
};

class BoundingBox{
  public:
    float min_x, max_x;
    float min_y, max_y;
    float min_z, max_z; 
    
    BoundingBox(); 
    BoundingBox(float, float, float, float, float, float); 
    
    int getLongestAxis(); 
    float getMidPoint(int); 
    bool intersect(Ray); 
};


class Shape{
  public:   
    BRDF brdf;
    virtual bool getIntersect(Ray, float*) = 0; 
    virtual bool getIntersect(Ray) = 0; 
    virtual Vector getNormal(Point) = 0; 
    virtual BoundingBox getBB() = 0; 
    
    bool isSphere;
    
};


class ShapeList{ //see AggregatePrimitive on website 
    public:
        vector<Shape*> allShapes; 
        bool checkIntersect(Ray, Point*, Shape*&, float);
        bool checkIntersect(Ray, float);
        
        BoundingBox getRootBox(); 
        
        ShapeList(); 
        ShapeList(vector<Shape*>); 
};


class Triangle: public Shape{
    public:
        vector<Point> vertices; 
        Point v1, v2, v3; 
        Matrix RtimesSinv;
    
        Triangle() {};
        Triangle(Point, Point, Point, Matrix); 
        
        bool getIntersect(Ray, float*); 
        bool getIntersect(Ray); 
        Vector getNormal(Point); //no need for Point b/c all normals are same on triangle
        BoundingBox getBB(); 

};

class AABB_Node{
  public:
        BoundingBox bb;
        AABB_Node* children[2];  
        vector<Triangle> containedTriangles; 
        
        AABB_Node(){}; 
        AABB_Node(vector<Triangle>*, int); 
        float CollisionTest(Ray, Triangle*, float*); 
        bool CollisionTest(Ray, float*); 
};


class Sphere: public Shape{
    public:
        Point center;
        float radius;
        
        Matrix R;
        Matrix Rinv;
        Matrix S;
        Matrix Sinv;
        Matrix invertTrans;
        
    
        Sphere() {}; 
        Sphere(Point, float); 
        Sphere(Point, float, Matrix, Matrix, Matrix, Matrix); 
        
        
  
        bool getIntersect(Ray, float*); 
        bool getIntersect(Ray); 
        Vector getNormal(Point); 
        BoundingBox getBB(); 
        
        Ray getLight(Ray);
        
      
};









