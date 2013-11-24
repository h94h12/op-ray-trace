#include "primitives.h"

//****************************************************
// brdf.h. get the color of a certain pixel with given lights
//****************************************************


//given a point on the unit sphere, return the normal vector 
Vector getNormal(Point p){
    return Vector(p.x, p.y, p.z).normalize();
}
//given normalized light and normal vectors, return normal reflection vector
Vector getReflection(Vector l, Vector n){
    l = l.normalize();
    n = n.normalize(); 
    double dp = l.dotProduct(n); 
    Vector v =  (n * 2 * dp) - l; 
    return v.normalize(); 
    
    
}
//given point on sphere and location of light, return normal direction vector to light
Vector getLight(Point p, Point l){
    Vector v = l.subtract(p); 
    return v.normalize();
    
}

