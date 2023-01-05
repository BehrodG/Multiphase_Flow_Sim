//
//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
//
#ifndef PFMedProp_Arrays 
#define PFMedProp_Arrays

// 
         //SAVE
// 
// ----------
// ...... Derived-type variables     
// ----------
// 

struct PFMedium
{
// 
	int   NumF_RP,     NumF_CP        ;// The number of the capillary function
// 
	float   DensG,  Poros,  SpcHt,  KThrW,   KThrD   ;// The thermal conductivity of the dry PF medium (W/m/C)
// 
	float   Compr,    Expan,   Tortu,    Klink   ;// The Kinkenberg parameter (1/Pa)
// 
	float   RGrain, RSpace,   SArea,      NVoid,   VVoid   ;// The volume of each interstitial void (m^3)
// 
	float DIMENSION(3)   Perm   ;// Intrinsic permeabilities along principal axes (m^2)
// 
	float DIMENSION(7)   ParRP  ;// Parameters of the relative permeability functions
// 
	float DIMENSION(7)   ParCP   ;// Parameters of the capillary pressure functions
// 
}
// 
// 
// 

struct PFMedium1
{
// 
	float   RGrain,   SArea,   RVoid,   NVoid,   VVoid   ;// The volume of each void (m^3)
// 

};
// 
// ----------
// ...... Derived-type arrays     
// ----------
// 

PFMedium    media;
// 
PFMedium1   PoMed;
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 
float XIN,YIN;
// 
// ----------
// ...... Double precision arrays     
// ----------
// 
float RPD,CPD
// 
// ----------
// ...... int allocatable arrays     
// ----------
// 
int StateI;
int StateIndex;
// 
// ----------
// ...... char allocatable arrays     
// ----------
// 
char  MAT[5],Rk_name[5]
// 
// 

#endif

