#ifndef BONDCOLLECTION_H
#define BONDCOLLECTION_H

#include "classdec.h"
#include "bonds.h"


class Bondcollection  //Will store bonds for whole system as a vector<Bondcollection>.  Each platelet has a Bondcollection.
{
	public:
	Bondcollection();
    ~Bondcollection();
     
    vector<Bonds> bonds;      //list of bonds on a platelet
    int weakbondcount;        //count of active bonds for platelet
	int strongbondcount;

    void ResetBondcollection();
    bool CheckBondCountsonPlatelet();
    void SetBondCountsonPlatelet();
    void RemoveBond(int );
};



#endif
