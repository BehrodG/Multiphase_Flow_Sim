//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
////
#ifndef Variable_Arrays 
#define Variable_Arrays

// 
         //SAVE
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 

float    X,DX,DELX;

float    R,DOLD;
// 
// ----------
// ...... Double precision arrays     
// ----------
// 

float DEP[10];   //float DIMENSION(10)   DEP //;

float XX[10];//float DIMENSION(10)   XX//;
// 
// ----------
// ...... Derived-type variables     
// ----------
// 

struct SecParam
{
//            float DIMENSION(:),   POINTER   p_Satr // Phase saturation
//            float DIMENSION(:),   POINTER   p_KRel // Phase relative permeability
//            float DIMENSION(:),   POINTER   p_Visc // Phase viscosity (Pa.s)
//            float DIMENSION(:),   POINTER   p_Dens // Phase density (kg/m^3)
//            float DIMENSION(:),   POINTER   p_Enth // Phase specific enthalpy (J/Kg) 
//            float DIMENSION(:),   POINTER   p_PCap // Phase capillary pressure (Pa)
//            float DIMENSION(:,:), POINTER   p_MasF // Mass fraction of components in the phases
//            float DIMENSION(:),   POINTER   p_DiFA // Coefficient 1 for multicomponent bin. diffusion
//            float DIMENSION(:),   POINTER   p_DiFB // Coefficient 2 for multicomponent bin. diffusion
//            float DIMENSION(:),   POINTER   p_AddD // Additional data 
//            float                            TemC   // Temperature (C)
	float *p_Satr; // Phase saturation
	float *p_KRel; // Phase relative permeability
	float *p_Visc; // Phase viscosity (Pa.s)
	float *p_Dens; // Phase density (kg/m^3)
	float *p_Enth; // Phase specific enthalpy (J/Kg) 
	float *p_PCap; // Phase capillary pressure (Pa)
	float *p_MasF; // Mass fraction of components in the phases
	float *p_DiFA; // Coefficient 1 for multicomponent bin. diffusion
	float *p_DiFB; // Coefficient 2 for multicomponent bin. diffusion
	float *p_AddD; // Additional data 
	float  TemC;   // Temperature (C)
};
// 
// ----------
// ...... Derived-type arrays     
// ----------
// 

SecParam *Cell_V; //SecParam ALLOCATABLE, DIMENSION(:,:)   Cell_V;
// 
// 
#endif

