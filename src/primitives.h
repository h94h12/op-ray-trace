//****************************************************
// primitives.h
// Contains:
// - Point
// - Ray
// - Vector
// - Light
// - Color
// - Camera 
// - MatrixStack
// - Matrix 
//****************************************************

#ifndef LIGHTS_H
#define LIGHTS_H

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


//****************************************************
// CLASS DECLARATIONS AND CONSTRUCTORS 
//****************************************************
class Vector; 
class Point{ 
public:
    Point() {};
    Point(float, float, float); 

    float x, y, z;
    Vector subtract(Point); 
    Point operator + (Vector); //position + direction = position 
    Point operator + (Point); // point addition; shortcut for camera's ray formula
    Point operator * (float); 
    Vector operator - (Point);  
    
    
    
};

    
    
class Vector{
 public:
    Vector() {};  
    Vector(float, float, float);
    Vector(Point); 
    Vector(Point, Point); 
    
    float dx, dy, dz;
    float mag;
    
    Vector normalize();
    float dotProduct(Vector);
    Vector negative();
    Vector crossProduct(Vector);
    
    //allow for scalar * and vector +
    Vector operator * (float);
    Vector operator + (Vector); 
    Vector operator - (Vector); 
    
    
};




class Ray{
public:
    Ray() {}; 
    Ray(float, float, float, float, float, float);
    Ray(Point, Point); 
    Ray(Point, Vector);

    Vector direction; 
    Point origin; 
   
    
};


/*
 *  Color class keeps track of r g b values, which are [0, 1]
 */ 
 class Color{
     public:
        Color() {};
        Color(float, float, float);
        
        float r, g, b;
        
        Color operator + (Color);  
        Color operator * (Color); 
        Color operator * (float); 
        void operator += (Color); 
        
        Color clone();
 };

/*
 *  Light class describes a light source that extends in all directions with rgb value from point x y z
 */ 

class Light {
  public:
    Light() {}; 
    Light (float, float, float, float, float, float, bool); 

    void initPos(float, float, float);
    void initRGB(float, float, float); 
    Color color; 
    Point source; 
    float x, y, z, r, g, b;
	bool directional;

};

class Matrix{
    public:
        vector<vector<float> > m; 
    Matrix();
    
    // Instantiator shortcut for translation, scaling, or rotating (depends upon input char)
    Matrix(char, float, float, float, float); 
    
    Matrix operator * (Matrix); 
    Point operator * (Point);
    Vector operator * (Vector);
    Ray operator * (Ray);
    
    Matrix invert();
    void debug();
    Matrix clone();
    
    Vector vectorTimesM(Vector); 
};

class MatrixStack{
    public:
    MatrixStack(); 
    
    vector<Matrix> stack;
    vector<Matrix> stackT;
    vector<Matrix> stackS;
    vector<Matrix> stackR;
    Matrix product;
    Matrix productT; // product of all the T matrices
    Matrix productS;
    Matrix productR;
    
    void push(); 
    void pop();
    void addTransform(Matrix, char);
     
};







 
#endif
