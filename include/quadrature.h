#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "classdec.h"
#include "tri.h"

// 2D quadrature
struct Quad2D {
    int _n;
    double *_w, *_x, *_y;

    int n() const { return _n; }
    double w(int i) const { return _w[i]; };
    double x(int i) const { return _x[i]; };
    double y(int i) const { return _y[i]; };

    double& w(int i) { return _w[i]; };
    double& x(int i) { return _x[i]; };
    double& y(int i) { return _y[i]; };

    Quad2D() 
    {
        _n = 0;
	_w = _x = _y = NULL;
    }
    
    ~Quad2D() 
    {
        deallocateMemory();
    }

    void allocateMemory(int n) 
    {
        _n = n;
        _w = new double[n];
        _x = new double[n];
        _y = new double[n];
    }

    void deallocateMemory() 
    {
        _n = 0;
	if (_w) delete [] _w;
	if (_x) delete [] _x;
	if (_y) delete [] _y;
    }
};


namespace quadrature 
{
    extern int _inited;
    extern Quad2D _Q_TRI_1;
    extern Quad2D _Q_TRI_3;
    extern Quad2D _Q_TRI_7;
    extern Quad2D _Q_QUAD_3;
    extern Quad2D _Q_QUAD_4;

    void _init();
    Quad2D& select_rule_2d(const char*);
};


int gauleg(double x1, double x2, int n, double *x, double *w);
void calcDuffyQuadPoint(int n, 
	vector<double> &sq, vector<double> &tq, vector<double> &wq);
void calcDuffyQuadPoint(double s0, double t0, int n, 
	vector<double> &sq, vector<double> &tq, vector<double> &wq);

void calcLogQuadPoint(double dist, double h, double thmin, double thmax,
	int n, vector<double> &sq, vector<double> &tq, vector<double> &wq);
void calcLogQuadPoint(const double *x0, const Tri &T,
	int n, vector<double> &sq, vector<double> &tq, vector<double> &wq);

#endif
