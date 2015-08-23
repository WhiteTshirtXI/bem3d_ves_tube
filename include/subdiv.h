#ifndef SUBDIV_H
#define SUBDIV_H

#include "classdec.h"
#include "mathfunc.h"

namespace subdiv {
    void divideMesh(const MArray<double,2> &x, const MArray<int,2> &f2v,
    		MArray<double,2> &xnew, MArray<int,2> &f2v_new);

    void calcFunc_666(double v, double w, 
	    MArray<double,1> &f, MArray<double,2> &df, MArray<double,3> &ddf);

    void calcFunc_Recursive(int N0, int N1, int N2, double v, double w, 
	    MArray<double,1> &f, MArray<double,2> &df, MArray<double,3> &ddf,
	    int level=0, int min_level=-10000);

    void calcFunc(int N0, int N1, int N2, double v, double w, 
	    MArray<double,1> &f, MArray<double,2> &df, MArray<double,3> &ddf, 
	    bool updateDb=false);

    void buildFuncDataBase(const Quad2D &);

    void buildControlPoints(int nvert, int nface, const int (*f2v)[3],
    		vector<vector<int> > &ctrlPoints);

    //**********************************************************************
    // Key to the Loop subdiv function database
    struct DbKey {
        int _N0, _N1, _N2;
	double _v, _w;

	DbKey(int N0, int N1, int N2, double v, double w) :
		_N0(N0), _N1(N1), _N2(N2), _v(v), _w(w) { }
    };

    // Values of the Loop subdiv functions
    typedef vector<double> DbVal;

    struct DbComp {
	bool operator() (const DbKey& key0, const DbKey& key1) const
	{ 
	    if (key0._N0 == key1._N0 &&
		key0._N1 == key1._N1 &&
		key0._N2 == key1._N2 &&
		fabs(key0._v - key1._v) < 1.E-10 &&
		fabs(key0._w - key1._w) < 1.E-10) return false;

	    int v0 = int(9999*key0._N0 + 888*key0._N1 + 77*key0._N2 + 10000*key0._v + 100*key0._w);
	    int v1 = int(9999*key1._N0 + 888*key1._N1 + 77*key1._N2 + 10000*key1._v + 100*key1._w);
	    return (v0 < v1);
        }
    };

    typedef map<DbKey,DbVal,DbComp> Db;
};

#endif
