#include "cxxheaders.h"
#include "point.h"

Point::Point() 
{
    x[0] = x[1] = x[2] = 0.0;
    pp = NULL;
    mesh = NULL;
    indx = -1;
    Gindx = -1;
}


Point::Point(double _x, double _y, double _z) 
{
    x[0] = _x;    
    x[1] = _y;    
    x[2] = _z;
}


Point::Point(const Point &pt) 
{
    x[0] = pt.x[0];
    x[1] = pt.x[1];
    x[2] = pt.x[2];
}


Point& Point::operator=(const Point &pt) 
{
    x[0] = pt.x[0];
    x[1] = pt.x[1];
    x[2] = pt.x[2];

    return *this;
}


Point::~Point() 
{ }
