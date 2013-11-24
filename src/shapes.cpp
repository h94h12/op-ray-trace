#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

#include <cstdlib>

#include <time.h>
#include <math.h>

#include "shapes.h"

BRDF::BRDF(){
    kd = Color(0.2, 0.2, 0.2); // diffuse term
    ks = Color(0, 0, 0); // specular term
    ke = Color(0.2, 0.2, 0.2); // ambient/glowing term
    kr = 0.3; // reflection
    sp = 20; 
}
BRDF::BRDF(Color k_d, Color k_s, Color k_e, float k_r, int s_p){
    kd = k_d, ks = k_s, ke = k_e, kr = k_r, sp = s_p; 
}
// a deep clone of the object
BRDF BRDF::clone(){
    BRDF n;
    n.kd = Color(kd.r, kd.g, kd.b);
    n.ks = Color(ks.r, ks.g, ks.b);
    n.ke = Color(ke.r, ke.g, ke.b);
    n.kr = kr;
    n.sp = sp;
    return n;
}

BoundingBox::BoundingBox(){
      min_x = 0, max_x = 0, min_y = 0, max_y = 0, min_z = 0, max_z = 0; 
}

BoundingBox::BoundingBox(float minx, float maxx, float miny, float maxy, float minz, float maxz){
    min_x = minx, max_x = maxx, min_y = miny, max_y = maxy, min_z = minz, max_z = maxz; 
}

int BoundingBox::getLongestAxis(){
    float xd = abs(max_x - min_x); 
    float yd = abs(max_y - min_y);
    float zd = abs(max_z - min_z);
    float m = max(xd, max(yd, zd));
    if(m == xd)
        return X_AXIS;
    else if (m == yd)
        return Y_AXIS;
    return Z_AXIS; 
    
}

float BoundingBox::getMidPoint(int axis){
    if(axis == X_AXIS){
        
        return (max_x + min_x)/2; 
    }
    if(axis == Y_AXIS){
        
        return (max_y + min_y)/2;  
    }
    else{
        
         return (max_z + min_z)/2; 
    }
    
}



Sphere::Sphere(Point c, float r){
    brdf = BRDF(); 
    center = c;
    radius = r; 
    //trans = Matrix(); // Identity by default
    //invertTrans = Matrix(); // Identity by default
    R = Matrix();
    Rinv = Matrix();
    S = Matrix();
    Sinv = Matrix();
    invertTrans = Matrix();
     
    isSphere = true; 
    
}
Sphere::Sphere(Point c, float r, Matrix rm, Matrix rinv, Matrix s, Matrix sinv){
    brdf = BRDF(); 
    center = c;
    radius = r; 
    R = rm;
    Rinv = rinv;
    S = s;
    Sinv = sinv;
    invertTrans =  Sinv * Rinv;
    //Matrix trans = R * S;

    
    
    isSphere = true; 
    
}

Vector Sphere::getNormal(Point p){
    //return (invertTrans * Vector(trans * center, trans * p)).normalize(); 
    //return (invertTrans * Vector(center, p)).normalize(); 
    
    
    return (R * Sinv * Rinv * Vector(center, p)).normalize(); 
    //return (Vector(center, p)).normalize(); 
    
    
    /*in order to the the surface normal in the coordinate system your elipsoid is in, 
     * take the vector from the sphere centre to the point of collision in transformed space,
     * and multiuply it with the adjunct (inverse transpose) elipsoid transformation matrix.
     *  then normalize if required. */
     
     //http://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
     //normal scales as INVERSE of object scaling, but rotates the same 
     
     
     
}   

BoundingBox Sphere::getBB(){
    return BoundingBox(center.x - radius,
                       center.x + radius,
                       center.y - radius,
                       center.y + radius,
                       center.z - radius,
                       center.z + radius); 
                       
    
}

    
bool Sphere::getIntersect(Ray rayinp, float* t){
	Ray ray = invertTrans * rayinp;

    Vector AC = ray.origin.subtract(invertTrans * center); 
    float temp1 = pow(AC.dotProduct(ray.direction), 2); 
    float temp2 = AC.dotProduct(AC) - pow(radius, 2); 
    float vdot = ray.direction.dotProduct(ray.direction);
    temp2 *= vdot;
	//
    if (temp1 - temp2 < 0)
        return false; 
    else {
		float temp3 = -2 * (AC.dotProduct(ray.direction)); 
		temp1 *= 4;
		temp2 *= 4;

		float soln = min(temp3 + sqrt(temp1 - temp2), temp3 - sqrt(temp1 - temp2)) / (2.0 * vdot);
		*t = soln; // soln is distance in object space 
		
		//Point po = ray.origin + (ray.direction * soln);
		//Point pw = /*trans * po;
		//*t = sqrt(pow(rayinp.origin.x - pw.x, 2) + pow(rayinp.origin.y - pw.y, 2)  + pow(rayinp.origin.z - pw.z, 2));
		
		//cout << soln << "__" << *t << endl;
		
		return true; 
    }
} 

// Deprecated
Ray Sphere::getLight(Ray L){
       // return (invertTrans * Vector(center, p)).normalize(); 
    //without regard to Transformations, light is Ray(bias, lightray,direction); 
    
    return Ray(L.origin,  L.direction);  
    //return Ray(L.origin, trans* L.direction);  
    
}



bool Sphere::getIntersect(Ray rayinp){
	Ray ray = invertTrans * rayinp;
	
	//ray.direction.normalize();
	
    Vector AC = ray.origin.subtract(invertTrans * center); 
    float temp1 = pow(AC.dotProduct(ray.direction), 2); 
    float temp2 = AC.dotProduct(AC) - pow(radius, 2); 
    float vdot = ray.direction.dotProduct(ray.direction);
    temp2 *= vdot;
    if (temp1 - temp2 < 0)
        return false; 
    return true; 
}


Triangle::Triangle(Point v_1, Point v_2, Point v_3, Matrix m){
    vertices.push_back(v_1);
    vertices.push_back(v_2);
    vertices.push_back(v_3);
    
    RtimesSinv = m;
    
    brdf = BRDF(); 
    
    v1 = v_1;
    v2 = v_2; 
    v3 = v_3; 
    
    isSphere = false; 
}

BoundingBox Triangle::getBB(){
    //don't return bounding box with 0 volume
    
   float xmin = min(min(v1.x, v2.x), v3.x);
   float xmax = max(max(v1.x, v2.x), v3.x);
   float ymin = min(min(v1.y, v2.y), v3.y);
   float ymax = max(max(v1.y, v2.y), v3.y);
   float zmin = min(min(v1.z, v2.z), v3.z);
   float zmax = max(max(v1.z, v2.z), v3.z);
   
    
    return BoundingBox(xmin, xmax,
                       ymin, ymax,
                       zmin, zmax); 
                       
}



Vector Triangle::getNormal(Point p){
    Vector v13 = v1 - v3; 
    Vector v23 = v2 - v3; 
    
    //return (RtimesSinv * v13.crossProduct(v23)).normalize();
    return v13.crossProduct(v23).normalize();
}


bool Triangle::getIntersect(Ray ray, float* t){ 
   Vector u = v1 - v3;
   Vector v = v2 - v3; 
   Vector N = u.crossProduct(v).normalize(); 

   if (N.dotProduct(ray.direction) == 0)
        return false; 
   float d = N.dotProduct(Vector(v1.x, v1.y, v1.z));  
   float soln = (d - (N.dotProduct(Vector(ray.origin.x, ray.origin.y, ray.origin.z)))) / N.dotProduct(ray.direction); 
   if (soln < 0 )
        return false; 
   Point interpt =  Point(ray.origin.x + soln * ray.direction.dx, 
                    ray.origin.y + soln * ray.direction.dy,
                    ray.origin.z + soln * ray.direction.dz);

   float surf3 = (v3 - v2).crossProduct(interpt - v2).dotProduct(N); 
   float surf2 = (v2 - v1).crossProduct(interpt - v1).dotProduct(N); 
   float surf1 = (v1 - v3).crossProduct(interpt - v3).dotProduct(N); 
   if(surf3 >= 0 && surf2 >= 0 && surf1 >= 0){ //intersection point inside all the edges of triangle
       *t = soln; 
       return true;
   }
   return false; 
}




bool Triangle::getIntersect(Ray ray){ //used algorithm from http://geomalgorithms.com/a06-_intersect-2.html

	float t; 
	return getIntersect(ray,  &t);    
}

ShapeList::ShapeList(){
    allShapes = vector<Shape *> ();  
}

ShapeList::ShapeList(vector<Shape*> shapes){
    allShapes = shapes; 
}

//find the closest shape the ray interects (FOR POSITIVE T), and at which point
//return false if no intersection
bool ShapeList::checkIntersect(Ray ray, Point* p, Shape*& sh, float tmax){
    float t = -1;
    float min = LARGE_NUM; 
    for(int i = 0; i < allShapes.size(); i++){
        allShapes[i]->getIntersect(ray, &t); 
        if (t < min && t > 0){
            min = t;
            sh = allShapes[i];  //buggy?
        }
    }
    if(min == LARGE_NUM || min > tmax) //no intersection
        return false;
    *p = Point(ray.origin.x + min * ray.direction.dx, 
				ray.origin.y + min * ray.direction.dy,
				ray.origin.z + min * ray.direction.dz); 

    return true; 
}
bool ShapeList::checkIntersect(Ray ray, float tmax){
    float t = -1;
    Shape* closest; 
    for(int i = 0; i < allShapes.size(); i++){
       allShapes[i]->getIntersect(ray, &t);
       if (t > 0 && t < tmax) //there is an intersection 
            return true; 
    }
    
    return false;
}

BoundingBox ShapeList::getRootBox(){
    float minx, miny, minz = LARGE_NUM; 
    float maxx, maxy, maxz = -LARGE_NUM; 
    BoundingBox bb; 
    
    for(int i = 0; i < allShapes.size(); i++){
        bb = allShapes[i]->getBB(); 
        minx = fmin(bb.min_x, minx);
        miny = fmin(bb.min_y, miny);
        minz = fmin(bb.min_z, minz);
        
        maxx = fmax(bb.max_x, maxx);
        maxy = fmax(bb.max_y, maxy);
        maxz = fmax(bb.max_z, maxz);
    }
    
    return BoundingBox(minx, maxx, miny, maxy, minz, maxz); 
    
}

bool BoundingBox::intersect(Ray r){ 
    if (!this) return false;  //how are you getting a null BB???
    
    
    float xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0; 
    float xd = r.direction.dx;
    if (xd == 0) xd = 0.00001;
    
    float yd = r.direction.dy;
    if (yd == 0) yd = 0.00001;
    
    float zd = r.direction.dz;
    if (zd == 0) zd = 0.00001;
    
    float a_x = 1 / xd;
    if(a_x >= 0){
        xmin = a_x * (min_x - r.origin.x);
        xmax = a_x * (max_x - r.origin.x);   //cannot access memory at max_x
    }
    else{
        xmin = a_x * (max_x - r.origin.x);
        xmax = a_x * (min_x - r.origin.x); 
        
    }
    
    float a_y = 1 / yd;
    if(a_y >= 0){
        ymin = a_y * (min_y - r.origin.y);
        ymax = a_y * (max_y - r.origin.y); 
    }
    else{
        ymin = a_y * (max_y - r.origin.y);
        ymax = a_y * (min_y - r.origin.y); 
        
    }
   
    if((xmin > ymax) || (ymin > xmax))
        return false;
    if(ymin > xmin)
        xmin = ymin;
    if(ymax < xmax)
        xmax = ymax;
      
    float a_z = 1 / zd;
    if(a_z >= 0){
        zmin = a_z * (min_z - r.origin.z);
        zmax = a_z * (max_z - r.origin.z); 
    }
    else{
        zmin = a_z * (max_z - r.origin.z);
        zmax = a_z * (min_z - r.origin.z); 
    }
    if((xmin > zmax) || (zmin > xmax))
        return false;
    if(zmin > xmin)
        xmin = zmin;
    if(zmax < xmax)
        xmax = zmax; 
    return (xmax > 0); 
   
   
}
BoundingBox getBB(vector<Triangle>* tris){
    float minx, miny, minz, maxx, maxy, maxz;  
    minx  = miny = minz = LARGE_NUM; 
    maxx = maxy = maxz = -LARGE_NUM; 
    BoundingBox bb = BoundingBox(); 
    
    for(int i = 0; i < tris->size(); i++){
        bb = tris->at(i).getBB(); 
        minx = fmin(bb.min_x, minx);
        miny = fmin(bb.min_y, miny);
        minz = fmin(bb.min_z, minz);
        
        maxx = fmax(bb.max_x, maxx);
        maxy = fmax(bb.max_y, maxy);
        maxz = fmax(bb.max_z, maxz);
    }
    
    return BoundingBox(minx, maxx, miny, maxy, minz, maxz); 
    
}
//simple method for intersection, just need to know if it intersects or not
bool checkTriIntersect(vector<Triangle>* tris, Ray ray, float* tmax){
   float t = -1;
    for(int i = 0; i < tris->size(); i++){
       tris->at(i).getIntersect(ray, &t);
       if (t > 0 && t < *tmax) //there is an intersection 
       {
            *tmax = t; 
            return true; 
        }
    }
    
    return false; 
}
//more detailed intersection test, also gives point of intersection
float checkTriIntersect(vector<Triangle>* tris, Ray ray, Triangle* tr, float* tmax){
    float t = -1;
    float min = *tmax; 
    for(int i = 0; i < tris->size(); i++){
        tris->at(i).getIntersect(ray, &t); 
        if (t > 0 && t < min){
            min = t;
            *tr = tris->at(i);  
        }
    }
    if(min == LARGE_NUM || min > *tmax) //no intersection
        return NO_INTERSECTION; 
    *tmax = min; 
    return min; 

    
}




//make tree by passing in allshapes to get root
//adapted from http://www.flipcode.com/archives/Dirtypunks_Column-Issue_05_AABB_Trees_Back_To_Playing_With_Blocks.shtml
AABB_Node::AABB_Node(vector<Triangle>* tris, int depth){ 

    bb = getBB(tris); 
     if(tris->size() >  3 && depth < 2){ 
        int axis = bb.getLongestAxis(); 
        vector<Triangle> left = vector<Triangle>(); 
        vector<Triangle> right = vector<Triangle>(); 
     
        float mid = bb.getMidPoint(axis); 
        //divide this bb by midpoint along longest axis, sort objects
        for(int i = 0; i < tris->size(); i++){
            if(tris->at(i).getBB().getMidPoint(axis) < mid)
                left.push_back(tris->at(i));
            else
                right.push_back(tris->at(i)); 
        }
        if(left.size() > 0){
           if(left.size() == tris->size()) children[0] = new AABB_Node(&left, depth + 1); 
           else children[0] = new AABB_Node(&left, depth);   
        } 
        else children[0] = 0; 
        if(right.size() > 0){
            if(right.size() == tris->size()) children[1] = new AABB_Node(&right, depth + 1); 
            else children[1] = new AABB_Node(&right, depth); 
        } 
        else children[1] = 0; 
        
    }
    else{ //make leaf node 
        children[0] = 0; 
        children[1] = 0; 
    }
    containedTriangles = *tris; 
    
}
//simpler test, don't need to know where ray intersects with shapes
bool AABB_Node::CollisionTest(Ray ray, float* limit){
     if(!bb.intersect(ray)){
        return false; 
    }

       if(children[0] || children[1]){
            if(children[0] && children[0]->CollisionTest(ray, limit))
                return true;
            if(children[1] && children[1]->CollisionTest(ray, limit))
                return true; 
            return false; 
        }
    
        
        //do intersection with shapes at this node 
    return checkTriIntersect(&containedTriangles, ray, limit); 
    
    
    
}
//return t value for ray, calculate point of intersection with t value
// if returns 0, no intersection
float AABB_Node::CollisionTest(Ray ray, Triangle* tri, float* limit){
    if(!bb.intersect(ray)){
        return NO_INTERSECTION; 
    }
    if(!children[0] && children[1]) return children[1]->CollisionTest(ray, tri, limit);
    if(!children[1] && children[0]) return children[0]->CollisionTest(ray, tri, limit); 
    if(children[0] && children[1])  return min(children[0]->CollisionTest(ray, tri, limit), children[1]->CollisionTest(ray, tri, limit)); 
        
    return checkTriIntersect(&containedTriangles, ray, tri, limit); 
}
