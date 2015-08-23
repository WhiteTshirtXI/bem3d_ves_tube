#ifndef BONDS_H
#define BONDS_H

#include "classdec.h"

//Remember to make these globals a namespace later
//Bond rate constants
const double bonds_k_on_weak = 5.;
const double bonds_k_off_weak = 5.;
const double bonds_k_on_strong = 5.;
const double bonds_k_off_strong = 0;
const double bonds_formationdist = 0.1;
const double bonds_springk = 15;  //Right now same k for weak/strong bonds
//const double bonds_springk_strong = 5;

//Maximum number of bonds to prevent platelet from pulling through wall and allow fixed array sizing
const int bondcollection_maxbonds = 8;



//Convention: the tail of the bond is on the platelet and the head is on the wall.

class Bonds  //Will store bonds for whole system as a vector<vector<Bonds>>.  Each rigid has a vector<Bonds>.
{
	public:
	Bonds();
    ~Bonds();
    
    //bool bondisactive;  //Active flag now merged to receptortype
	int bondreceptortype;  //Right now -17 = inactive, 0 = weak bond, 1 = strong bond.  
	               //This is an int not a bool for future use when there may be more types of bonds or drug interactions that affect the bonds
	int bondtailvert;  //Simplified version: bond lives at vertex of rigid.  Get rid of this later.               
	int bondtailtri;   //Triangle on platelet where bond lives
	double bondtailuv[2];  //Barycentric coordinate on platelet headtri where bond lives
	double bondheadx[3];   //Coordinate on wall where bond lives
};



#endif
