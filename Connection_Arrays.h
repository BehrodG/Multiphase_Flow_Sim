//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
///
#ifndef Connection_Arrays 
#define Connection_Arrays


// 
         //SAVE
// 
// ----------
// ...... Derived-type variables: Connection     
// ----------
// 

struct Connection
{
// 
	int   n1,n2;                                 // Numbers of elements in the connection
	int   ki;                                    // Number of the permeability index
// 

	char  nam1[8],nam2[8];                  // Element names in the connection
// 

	float   d1, d2;                          // Distance from interface (m)
	float   area;                           // Interface area (m^2)
	float   beta;                           // Cos(angle) with vertical
// 

	float   FluxH;                          // Heat flux (W)




// 
            //float DIMENSION(:), POINTER   FluxF   // Phase fluxes (kg/s)
            //float DIMENSION(:), POINTER   VelPh   // Phase velocities (m/s)


	float  *FluxF;   // Phase fluxes (kg/s)
	float  *VelPh;  // Phase velocities (m/s)
// 

};

// 
// ----------
// ...... Derived-type variables: VaporFlux     
// ----------
// 

struct VaporFlux
{
// 

	float  FluxVD;                // Diffusive vapor flux (W)
// 

	float  FluxVF[2];  // Vapor flux in the phases (kg/s)
// 

};
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 

float    sig;
// 
// ----------
// ...... int allocatable arrays     
// ----------
// 

int    NEX1_a,NEX2_a,IFL_con;
// 
// ----------
// ...... Derived-type arrays     
// ----------
// 

Connection  *conx;
// 

VaporFlux    *conV;
// 
#endif
