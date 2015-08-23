#include "cxxheaders.h"
#include "bondcollection.h"
#include "bonds.h"
#include "mathfunc.h"


Bondcollection::Bondcollection()
{ }


Bondcollection::~Bondcollection()
{ }


void Bondcollection::ResetBondcollection()
{//Initialize or Clear all bonds
	weakbondcount = 0;
	strongbondcount = 0;
	
	bonds.resize(bondcollection_maxbonds);
	for(int bonditer=0; bonditer < bondcollection_maxbonds; bonditer++)
	{
		//bonds[bonditer].bondisactive = false;
		bonds[bonditer].bondreceptortype = -17; //-17 is an uninitialized error default value
		bonds[bonditer].bondtailvert = -1;
		bonds[bonditer].bondtailtri = -2;
		bonds[bonditer].bondtailuv[0] = -3;
		bonds[bonditer].bondtailuv[1] = -4;
		bonds[bonditer].bondheadx[0] = -5;
		bonds[bonditer].bondheadx[1] = -6;
		bonds[bonditer].bondheadx[2] = -7;
		
	}
}

bool Bondcollection::CheckBondCountsonPlatelet() //Check count of weak and strong bonds to see if weakboundcount and strongbondcount are consistent with that's actually in the vector of bonds
{
	int observedweakbonds = 0;
	int observedstrongbonds = 0;
	for(int bonditer=0; bonditer < bondcollection_maxbonds; bonditer++)
	{
		if(bonds[bonditer].bondreceptortype !=-17)
		{
			if(bonds[bonditer].bondreceptortype == 0)
			{observedweakbonds++;}
			else if(bonds[bonditer].bondreceptortype == 1)
			{observedstrongbonds++;}
			//else if(bonds[bonditer].bondreceptortype == -17)
			//{cout << "BONDACTIVE17 ERROR: Bond " << bonditer << " is active but has no type!\n";}
			else //Fix this line if we add more types than weak/strong bonds:
			{cout << "BONDTYPEUNKNOWN ERROR: Bond " << bonditer << '-' << bonds[bonditer].bondreceptortype << " is neither weak nor strong nor initialized.\n";}
		}
	}
	if(weakbondcount != observedweakbonds)
	{cout << "WEAKBONDCOUNT ERROR: Bondcollection thinks there are " << weakbondcount << " weak bonds but there are really " << observedweakbonds << endl;}
	
	if(strongbondcount != observedstrongbonds)
	{cout << "STRONGBONDCOUNT ERROR: Bondcollection thinks there are " << strongbondcount << " strong bonds but there are really " << observedstrongbonds << endl;}
	
	return((weakbondcount == observedweakbonds) && (strongbondcount == observedstrongbonds));
}

void Bondcollection::SetBondCountsonPlatelet() //Instead of just checking, overwrite weakboundcount and strongbondcount with that's actually in the vector of bonds
{
	int observedweakbonds = 0;
	int observedstrongbonds = 0;
	for(int bonditer=0; bonditer < bondcollection_maxbonds; bonditer++)
	{
		if(bonds[bonditer].bondreceptortype !=-17)
		{
			if(bonds[bonditer].bondreceptortype == 0)
			{observedweakbonds++;}
			else if(bonds[bonditer].bondreceptortype == 1)
			{observedstrongbonds++;}
			//else if(bonds[bonditer].bondreceptortype == -17)
			//{cout << "SETBONDACTIVE17 ERROR: Bond " << bonditer << " is active but has no type!\n";}
			else //Fix this line if we add more types than weak/strong bonds:
			{cout << "SETBONDTYPEUNKNOWN ERROR: Bond " << bonditer << '-' << bonds[bonditer].bondreceptortype <<" is neither weak nor strong nor initialized.\n";}
		}
	}
	
	weakbondcount = observedweakbonds;
	strongbondcount = observedstrongbonds;
}



void Bondcollection::RemoveBond(int removeme)
{
	if(bonds[removeme].bondreceptortype==0)
	{weakbondcount=weakbondcount-1;}
	else if(bonds[removeme].bondreceptortype==1)
	{strongbondcount=strongbondcount-1;}
	else if(bonds[removeme].bondreceptortype==-17)
	{cout << "BONDALREADYINACTIVE ERROR: Bond" << removeme << " is already inactive!\n";}
	else
	{cout << "DELETEBONDTYPEUNKNOWN ERROR: Bond " << removeme << " is slated for removal but is neither weak nor strong.\n";}

	bonds.erase(bonds.begin()+removeme);
    bonds.resize(bondcollection_maxbonds);
    
	//bonds[bondcollection_maxbonds].bondisactive = false;
	bonds[bondcollection_maxbonds-1].bondreceptortype = -17; //-17 is an uninitialized error default value
	bonds[bondcollection_maxbonds-1].bondtailvert = -1;
	bonds[bondcollection_maxbonds-1].bondtailtri = -2;
	bonds[bondcollection_maxbonds-1].bondtailuv[0] = -3;
	bonds[bondcollection_maxbonds-1].bondtailuv[1] = -4;
	bonds[bondcollection_maxbonds-1].bondheadx[0] = -5;
	bonds[bondcollection_maxbonds-1].bondheadx[1] = -6;
	bonds[bondcollection_maxbonds-1].bondheadx[2] = -7;
}
