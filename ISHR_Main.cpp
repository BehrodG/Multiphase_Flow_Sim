// In the name of GOD

//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//*                                                                                            *
//*                                                                                            *
//*              Program of "Iranian Simulator for Hydrate Reservoirs"  Version 1.0            *
//*                                    28th March, 2014                                        *
//*                                    B. Gharedaghloo                                         *
//*                                                                                            *
//*                                                                                            *
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************




//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>                                                                     >
//>                                                                     >
//>     SHR_Main.cpp: Code unit including the main program              >
//>                      and all related routines                       >
//>                                                                     >
//>                                                                     >
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
//
//

#include <iostream>
#include <fstream>
#include <string>
#include <limits>
using namespace std;



// defintion for  Basic_Param

// 
         //SAVE
// 
// ----------
// ...... Double precision parameters
// ----------
// 

double   ZERO = 1.0e-6;
float   ZA   = 1.0e00;
// 
// ----------
// ...... Common int variables      
// ----------
// 
int   NEQ, NK, NPH, M_BinDif, M_add, Num_MobPhs;   
// 
int   No_CEN;
// 
int   NEQ1,NK1,NFLUX;
// 
int   MaxNum_Elem,NEL,NELA;
// 
int   MaxNum_Conx,NCON;
// 
int   MaxNum_PVar,Num_SecPar;
// 

int   NREDM,MNZ,NZ,LIRN,LICN,NCONUP,NCONDN,MNZP1;

int   NRWORK,NIWORK,lenw,leniw;
// 

int   MaxNum_SS,NOGN;
// 

int   MaxNum_Media,N_PFmedia;
// 

int   nactdi,nactd2;
// 

int   MaxNum_GasComp;  // Maximum number of gas components
// 
// ----------
// ...... char parameters
// ----------
// 

//char *GPhase_Com[6]; // Names of the gas components  

string *GPhase_Com; // Names of the gas components 
// 
// ----------
// ...... Common char variables
// ----------
// 

char FLG_con[6];

string EOS_Name;

char TITLE[80];
// 
// ----------
// ...... Common bool variables
// ----------
// 

bool   EX_FileVERS;
// 
// 
//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
///
//
//define Connection_Arrays


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

float    *sig;
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
//

//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
///
///
//define Diffusion_Param


// 
         //SAVE
// 
// ----------
// ...... Double precision arrays     
// ----------
// 
         
float    VBC;

float    FDIF;
// 
 
float   diffusivity[10][5];
// 
// ----------
// ...... Common double precision variables     
// ----------
// 

float   DIFF0,TEXP,BE;
// 
// ----------
// ...... Common int variables     
// ----------
// 

int   iddiag;
// 
// 
	  
//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
///


//define Solver_Param


//SAVE
// 
// ----------
// ...... Derived-type variables     
// ----------
// 

struct Element
{
// 
	int  MatNo  ; //int(KIND = 2)   MatNo  ;     // Element material #
	int  StateIndex ; //int(KIND = 1)   StateIndex ; // Element phase state index
	int  StatePoint;  //int(KIND = 1)   StatePoint;  // Element phase state point
// 
	char name[8]  ;     //char(LEN = 8)   name  ;     // Element name
// 
	float   vol  ;          // Element volume
	float   phi  ;          // Element porosity
// 
	float   P   ;           // Element pressure
	float   T   ;           // Element temperature
// 
};
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 

float    *AI,*AHT,*Pm;

float    *X_Coord;

float    *Y_Coord;

float    *Z_Coord;
// 
// ----------
// ...... Derived-type arrays     
// ----------
// 
 
Element  *elem;



//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//

//define EOS_Parameters
// 
//         //SAVE
// 
// ----------
// ...... int parameters
// ----------
// 
         
int MaxNum_MComp;   // Maximum number of mass components

int MaxNum_Equat;  // Maximum number of equations

int MaxNum_Phase;   // Maximum number of phases

int MaxNum_MobPh;   // Maximum number of mobile phases

int M_Add;



//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//


//define EOS_DefaultParam


//#include <string>
//#include <iostream>

//char EOS_state[3][11] = ( 'Gas', 'Aqu', 'AqG',  'IcG', 'AqH', 'IcH',  'AGH', 'AIG', 'AIH', 'IGH', 'QuP' );
//char EOS_state[3][11] = { {"Gas"}, {"Aqu"}, {"AqG"},  {"IcG"}, {"AqH"}, {"IcH"},  {"AGH"}, {"AIG"}, {"AIH"}, {"IGH"}, {"QuP"} );
	
string EOS_state[11]={"Gas", "Aqu", "AqG",  "IcG", "AqH", "IcH",  "AGH", "AIG", "AIH", "IGH", "QuP" };


//          CONTAINS


//void EOS_DefaultNum (string EOS_Name, int MaxNum_MComp, int MaxNum_Equat, int MaxNum_Phase, int MaxNum_MobPh, int M_Add);
void EOS_DefaultNum (string EOS_Name);







//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
//
//define GenControl_Param
// 
         //SAVE
// 
// ----------
// ...... Double precision arrays     
// ----------
// 

float DLT[100], TIS[100];
// 
// ----------
// ...... int arrays     
// ----------
// 

int   MOP[24];
// 
// ----------
// ...... Common double precision variables     
// ----------
// 

float   U_p;

float   TZERO,T_mat_solv;

float   TIMP1,TIMIN,TIMAX,REDLT,SUMTIM;

float   TSTART,DELAF,DELTMX,TIMOUT,GF;

float   DELT,DELTEN,DELTEX,FOR,FORD;

float   RE1,RE2,RERM,DFAC;

float   SCALE,WUP,WNR;
// 
// ----------
// ...... Common int variables     
// ----------
// 

int   NDLT,ITI,ITPR,ITCO;

int   KC,KCYC,ITER,ITERC;

int   INUM,IPRint, MCYC,MCYPR,MSEC;

int   IRPD,ICPD;

int   GenStIndx;
// 

int   KDATA,NER,KER,NOITE;

int   IGOOD,NST,IHLV;

int   MM;

int   IHALVE;
// 
// ----------
// ...... Common char variables     
// ----------
// 

//char  PermModif[4];
string PermModif;
char  ELST[8];

// 
// ----------
// ...... Common logical variables     
// ----------
// 

bool   CoordNeed,WellTabData,NoFloSimul,NoVersion;

bool   PrintOut_Flag,Converge_Flag,RandmOrdrIC;

bool   RadiHeat_Flag,CondHeat_Flag;


//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//


//define Input_Comments

// 
         //SAVE
// 
// ----------
// ...... Double precision arrays     
// ----------
// 

char *comm;

//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//



//define PFMedProp_Arrays

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
	float Perm[3]   ;// Intrinsic permeabilities along principal axes (m^2)
// 
	float ParRP[7]  ;// Parameters of the relative permeability functions
// 
	float ParCP[7]   ;// Parameters of the capillary pressure functions
// 
};
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

PFMedium    *media;
// 
PFMedium1   *PoMed;
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 
float *XIN,*YIN;
// 
// ----------
// ...... Double precision arrays     
// ----------
// 
float RPD[7],CPD[7];
// 
// ----------
// ...... int allocatable arrays     
// ----------
// 
int *StateI;
int *StateIndex;
// 
// ----------
// ...... char allocatable arrays     
// ----------
// 
//char  *MAT
string *MAT;
char *Rk_name;
// 
// 






//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************

//define Q_Arrays

// 
         //SAVE
// 
// ----------
// ...... Derived-type variables     
// ----------
// 

struct SourceSink
{
//

	char   name[5];
 
	char   name_el[8];
 
	char   type[5],dlv_file_name[5];
//
 
	int    el_num;
 
	int    n_TableP;
   
	int    typ_indx;
//
   
	float   rate_m,rate_v,rate_res,enth;
  
	float   pi,bhp,z_layer,grad;
//
   
	float  *frac_flo;
    
	float  *rate_phs;
//
	
};
// 
// 
// 
    
struct Tabular
{
         
	int   N_points;
          
	float  *p_TList, *p_QList;
         
	float  *p_HList;
		
};
// 
// 
// 
       
struct Tab_BHP
{
 
	int    idTab,NumG,NumH;

	//float DIMENSION(:), POINTER   p_bhp;

	float *p_bhp;

};
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 

float    *F1,*F2,*F3;
// 
// ----------
// ...... Derived-type arrays     
// ----------
// 

SourceSink  *SS;
// 
 
Tabular     *Well_TData;
// 
 
Tab_BHP     *WDel;





//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************


//define SELEC_Param


         //SAVE
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 
         float    FE;
// 
// ----------
// ...... Common int arrays    
// ----------
// 
         int   IE[16];
// 
// 





//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************



//define SolMatrix_Arrays

// 
         //SAVE
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 
float    CO;
float    *WKAREA;
float    RWORK;
// 
//float ALLOCATABLE, DIMENSION(:,:)   AB
float    *AB;
// 
// ----------
// ...... int allocatable arrays     
// ----------
// 

int    IRN,ICN;
int    JVECT;
// 

int    IIJJKK,NO;
int    IJKMM,IJKPP;
int    ISUMMM,ISUMPP;
// 

int  NX0;
int  NA0,NB0,NC0;
// 
// 





//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************





//define Solver_Param


// 
         //SAVE
// 
// ----------
// ...... Double precision parameters
// ----------
// 
         
double   seed = 1.0e-25;
// 
// ----------
// ...... Double precision variables
// ----------
// 
 
float   ritmax,closur;
// 
// ----------
// ...... int parameters
// ----------
// 
 
int iunit = 0, nvectr = 30;
// 
// ----------
// ...... Common int variables     
// ----------
// 

int   matslv,nmaxit,nnvvcc,iiuunn;
// 
// ----------
// ...... Common char variables     
// ----------
// 
 
char oprocs[2],zprocs[2];
 
char coord[4];
// 
// ----------
// ...... Common logical variables     
// ----------
// 
 
bool   SolverBlok; 

//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//


//#define TempStoFlow_Arrays

// 
         //SAVE
// 
// ----------
// ...... Double precision arrays     
// ----------
// 

float DEPU[10];
// 

float D[11][12];
// 

float F[11][23];
// 
// 



//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//


//define TimeSeriesP_Param

// 
         //SAVE
// 
// ----------
// ...... Derived types    
// ----------
// 

struct Observation_Elem
{
// 

	int  num[100];   // Element # of observation elements 
 
	char name[8][100];  // Names of observation elements
// 

};
// 
// 
// 

struct Observation_Conx
{
 
	int  num[100] ; // Connection # of observation connections
 
//	char name[16][100] ; // Names of observation connections
		
	string name[100] ; // Names of observation connections

};


// 
// 
// 

struct Observation_SS
{

	int    num[100];   // Source/sink # of observation sources/sinks
 
	char   name[5][100];  // Names of observation sources/sinks 

};

// 
// ----------
// ...... Derived-type variables     
// ----------
// 
 
Observation_Elem   obs_elem ;  // The observation element list

Observation_Conx   obs_conx ;  // The observation connection list
 
Observation_SS     obs_SS ;    // The observation SS lis
// 
// ----------
// ...... int variables
// ----------
// 

int   N_obs_elem,  N_obs_conx,      N_obs_SS ;      // The # of sources/sinks in the observation SS list
// 
// ----------
// ...... Logical variables
// ----------
// 
 
bool   obs_elem_Flag,   obs_conx_Flag,   obs_SS_Flag  ;   // Flag indicating that the observation SS list is not empty
// 



//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************



//#define Universal_Param

// ----------
// ...... Double precision parameters
// ----------
// 
         double GRPI  = 3.14159265358979324e0;
         double   R_gas = 8.31456e3;
// 
         double   Stef_Boltz = 5.6687e-8;
// 
// 



		 //***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
////

//define Variable_Arrays

// 
         //SAVE
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 

float    *X,*DX,*DELX;

float    *R,*DOLD;
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


// Define Hydrate_param 

// ----------
// ...... Double precision parameters
// ----------
// 
         float MW_CH4 = 16.04300;
// 
         float Tc_quad_hyd = 0.01;  // Quadruple point temperature
// 
// ----------
// ...... Double precision variables
// ----------
// 
         float P_quad_hyd;  // Quadruple point pressure
// 
         float MW_Hydrate;  // Molecular weight of the composite hydrate
// 
         float Initial_HydMass;  // Initial hydrate mass
// 
         float Area_Factor;      // Initial hydrate mass
// 
         float T_MaxOff,   C_MaxOff,   MW_Inhib,  D_Inhib,   H_InhSol,  DifCo_Inh;  // Diffusion coefficient of inhibitor (m^2/s)
// 
         float Activ_En,   InRate_k0;   // At input: intrinsic dissoc. rate constant [mol/(m^2.Pa.s)]
//                                      //    Later converted to [Kg/(m^2.Pa.s)] in "READ_HydrateData"
// 
// ----------
// ...... Derived-type variables     
// ----------
// 
         struct HydrProperties
		 {
			 int N_ThC;           // # of coefficients in the k_theta polynomial
			 int N_SpH;           // # of coefficients in the sp. heat polynomial
			 int N_Rho;           // # of coefficients in the density polynomial
// 
            float *p_ThC;           // Coefficients of the hydrate thermal conductivity 
//                                      // k_theta polynomial function [W/m/C^n, n=1,...,N_ThC]
            float *p_SpH;           // Coefficients of the hydrate specific 
//                                      // k_theta polynomial function [W/m/C^n, n=1,...,N_SpH]
            float *p_Rho;           // Coefficients of the hydrate density 
//                                      // polynomial function [kg/m^3/C^(n-1), n=1,...,N_Rho]
//                                      // Later converted to Kg/... in "READ_MainInputF"
		 };
		 

          struct HydrComposition{
			
			  int     NCom;   // The number of components of the ... 
//                                                      // ... composite hydrate
//			  char   nameG[6][2];  // The name of the hydrate-forming gases ...
//                                                      // ... in the composite hydrate
  			  string   nameG[2];  // The name of the hydrate-forming gases ...
//                                                      // ... in the composite hydrate

            
			  float hydrN[2], moleF[2], MolWt[2],   GasMF[2],  H2OMF[2];  // Mass fraction of H2O in the hydrate
// 
		  
			  
		  };

// 
// 
         HydrProperties GH;
         HydrComposition HCom;
		 

// 
// ----------
// ...... Integer parameters
// ----------
// 
         int F_EqOption;    // Flag for selection of hydrate Peq - Teq equilibrium function 
                                  // = 0: Moridis [2003]; = 1: Kamath (1984)
// 
// ----------
// ...... Logical parameters
// ----------
// 
         bool F_Kinetic;    // Flag indicating kinetic dissociation/formation
// 
// ----------
// ...... Double precision allocatable arrays     
// ----------
// 
         float *SArea_Hyd; // Area of hydration reaction
         float *Inh_WMasF;  // Mass fraction of inhibitor in injected stream
// 

//
//
//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***************************

     // MODULE Gas_DerivedType
 
//         SAVE
//! 
//! ----------
//! ...... Derived-type variables     
//! ----------
//! 
         struct GasParV{
			 string Name;
			 float AtomDV;
            
			 float TCrit,PCrit,VCrit,ZCrit;
            
			 float Omega,MolWt,DMom;
            
			 float A0,A1,A2,A3,A4;
             float Psi,A,B,T_B;
            

			 double LogK[5];  // Log10(K) regression coefficients
            
			 float MVol[6];   // Molar volume regression coefficients
            
			 float SOCo[6];  // Salting-out coefficient
		 } ;

//
         struct GasPar2
		 {
			 float TCrit,PCrit,ZCrit,VCrit;
           
			 float Omega,MolWt,DMom;
		 } ;
		 

// 
  




//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//      MODULE RefRGas_Param
 
//         USE Gas_DerivedType
 
//         SAVE
 
// ----------
//...... Double precision parameters
//----------
//
          const float Rgas = 8.314e3;

// ----------
//...... Integer parameters
// ----------
// 
          int NCom = 9;

// ----------
// ...... Derived-type parameters   
//

//...... NOTE: The Psi, A, B in the parameters
// ......    For CO2 : From Diamond and Atkinfiev (2003), Fluid Phase Equilib (in press) 
// ......    For rest: Table 1 and 2, Atkinfiev and Diamond, 2003, GCA (see above)
// 
// ...... NOTE: Salting-out coefficients
// ......    Langmuir (1997)  Aqueous Environmental Geochemistry, Prentice Hall, Table 4.5, p.144 
// ......    Modified because of s problem with the log form in the reference
// -------------
// 
// ----------
// 
//        TYPE(GasParV), PARAMETER :: Methane = GasParV
//     &       ("CH4   ", 2.514d1,                                     // Name, AtomDV
//     &         1.9056d+02,  4.600155d+06,  9.9000d-02,  2.880d-01,   // TCrit,PCrit,VCrit,ZCrit
//     &         1.1000d-02,  1.604300d-02,  0.0000d+00,  4.568d+00,   // Omega,MolWt,DMom,A0  
//     &        -8.9750d-03,  3.631000d-05, -3.4070d-08,  1.091d-11,   // A1,A2,A3,A4 
//     &        -0.1131d+00, -1.184620d+01,  1.48615d01,  111.6d0,     // Psi,A,B,T_B
//     &      (/-2.5764d+00, -1.3019d-02,    1.0924d-04, -2.9008d-07,  // Log10(K):  0-350 C
//     &         3.2774d-10 /),  
//     &      (/ 3.7076d+01,  2.1313d-03,   -2.8752d-05,  1.1936d-06,  // MolVol: 0-250 C 200 bar: SUPCRT92 + modified SLOP98
//     &        -7.3216d-09,  1.8923d-11 /),
//     &      (/ 1.64818d-1, -1.40166d-3,    1.32360d-5, -4.85733d-8,  // Salting-out coefficient: From Cramer
//     &        7.87967d-11, -5.52586d-14 /) )


		  GasParV Methane, Ethane, Propane, H2Sulfide, CarbonO2, Nitrogen, Oxygen, Water, Ethanol ;

/*		  // ["CH4   ", 2.514e1, 1.9056e+02,  4.600155e+06,  9.9000e-02,  2.880e-01, 1.1000e-02,  1.604300e-02,  0.0000e+00,  4.568e+00, -8.9750e-03,  3.631000e-05, -3.4070e-08,  1.091e-11, -0.1131e+00, -1.184620e+01,  1.48615e01,  111.6e0, (/-2.5764e+00, -1.3019e-02,    1.0924e-04,-2.9008e-07,3.2774e-10 /), (/ 3.7076e+01,  2.1313e-03,   -2.8752e-05,  1.1936e-06,  -7.3216e-09,  1.8923e-11 /),  (/ 1.64818e-1, -1.40166e-3,    1.32360e-5, -4.85733e-8, 7.87967e-11, -5.52586e-14 /) ];

		  //Ethane.Name="C2H6";
// 
//    Alternative for CH4:  25-300 C 350 bar: Regression from Hnedkovsky, MVOL = ...
//    (/ 3.5821E+01, 5.8603E-02,-1.3091E-03, 1.9022E-05,-8.6735E-08, 1.4126E-10 /) 
// 
         GasParV Ethane = ["C2H6  ", 4.566e1,3.0532e+02,  4.883865e+06,  1.4800e-01,  2.850e-01,9.9000e-02,  3.007000e-02,  0.0000e+00,  4.178e+00, -4.4270e-03,  5.660000e-05, -6.6510e-08,  2.487e-11,-0.6091e+00, -1.634820e+01,  2.00628e01,  184.6e0,(/-2.1993e+00, -1.8015e-02,    1.4569e-04, -4.0442e-07, 4.6619e-10 /),(/ 4.6452e+01,  2.3782e-01,   -2.8260e-03,  2.0155e-05, -7.1926e-08,  1.0645e-10 /),(/ 1.64818e-1, -1.40166e-3,    1.32360e-5, -4.85733e-8,7.87967e-11, -5.52586e-14 /) ];
		 GasParV Propane = ["C3H8  ", 6.618e1, 3.6983e+02,  4.2455175e+06,  2.0300e-01,  2.810e-01,   1.5200e-01,  4.4097000e-02,  0.0000e+00,  3.487e+00,  5.1310e-03,  6.0110000e-05, -7.8930e-08,  3.079e-11,  -1.1471e+00, -2.5387900e+01,  2.82616e01,  231.1e0,(/-2.3150e+00, -2.0996e-02,     1.7215e-04, -4.8268e-07,2.0620e-10 /), (/ 6.1525e+01,  2.8840e-01,    -3.4270e-03,  2.4452e-05, -8.7292e-08,  1.2927e-10 /), (/ 0.00e+00,    0.00e+00,       0.00e+00,    0.000e+00,   0.00e+00,    0.00e+00 /) ];
		 GasParV H2Sulfide = ["H2S   ", 2.752e1,  3.7340e+02,  8.9368650e+06,  9.8500e-02,  2.840e-01,  9.0000e-02,  3.4082000e-02,  9.0000e-01,  3.259e+00,   -3.4380e-03,  1.3190000e-05, -1.3310e-08,  0.488e-11, -0.2102e+00, -1.1230300e+01,  1.26104e01,  0.0e0,(/-6.7372e-01, -1.4245e-02,     7.9816e-05, -1.9861e-07, 5.5931e-10 /),   (/ 3.3065e+01,  8.3733e-02,    -9.9447e-04,  6.9408e-06,  -2.4299e-08,  3.4721e-11 /),(/ 2.9050e-1,  -1.57400e-4,    -4.62000e+1,  5.70500e-1,  -1.7770e-3,   0.00d+00   /) ];

// 
//    Alternative for H2S:  25-350 C 350 bar: Regression from Hnedkovsky, MVOL = ...
//    (/ 3.3521E+01, 6.8759E-02,-7.8118E-04, 9.2859E-06,-3.8445E-08, 6.1362E-11 /) 
// 
         GasParV CarbonO2 = ["CO2   ", 2.690e1, 3.0412e+02,  7.3764600e+06,  9.4000e-02,  2.740e-01, 2.2500e-01,  4.4010000e-02,  0.0000e-00,  3.539e+00, 1.3560e-03,  1.5020000e-05, -2.3740e-08,  1.056e-11, -0.0880e+00, -9.3134000e+00,  1.15477e01,  0.0e0, (/-1.1161e+00, -1.6023e-02,     9.5730e-05, -2.3733e-07, 2.3441e-10 /), (/ 2.9844e+01,  1.3722e-01,    -1.6281e-03,  1.0898e-05, -3.6679e-08,  4.8444e-11 /), (/-1.0312e00,   1.28060e-3,     2.55900e+2,  4.44500e-1, -1.6060e-3,   0.00e+00   /) ];

// 
//    Alternative 1 for CO2:  25-350 C 350 bar: Regression from Hnedkovsky, MVOL = ...
//    (/ 3.3101E+01, 6.1728E-03, 3.2205E-06, 6.3807E-06,-3.6252E-08, 6.8151E-11 /)
// 
//    Alternative 2 for CO2:  25-350 C 350 bar: Regression of J.E. Garcia, LBNL, 2001, MVOL = ...
//    (/ 3.7510E+01,-9.5850E-02,+8.7400E-04,-5.0440E-07,-4.4956E-19, 6.3527E-22 /) 

// 
         GasParV Nitrogen = ["N2    ", 1.850e1,        1.2620e+02,  3.3943875e+06,  8.9500e-02,  2.900e-01,   3.7000e-02,  2.8014000e-02,  0.0000e-00,  3.539e+00,    -0.2610e-03,  0.0070000e-05,  0.1570e-08, -0.099e-11,   -0.0320e+00, -1.1538000e+01,  1.46278e01,  0.0e0,  (/-2.9831e+00, -9.9098e-03,     8.5586e-05, -2.4931e-07,   3.0552e-10 /),    (/ 3.0349e+01,  1.3953e-01,    -1.6592e-03,  1.2163e-05,    -4.4426e-08,  6.8446e-11 /),   (/ 8.68589e-2,  0.00e+00,       0.00e+00,    0.000e+00,     0.00e+00,    0.00e+00 /) ]
// 
         GasParV Oxygen =   ["O2    ", 1.630e1,     1.5458e+02,  5.0459850e+06,  7.3400e-02,  2.880e-01,    2.1000e-02,  3.1999000e-02,  0.0000e-00,  3.630e+00,    -1.7940e-03,  0.6580000e-05, -0.6010e-08,  0.179e-11,     0.0260e+00, -9.7540000e+00,  0.00000e00,  0.0e0,  (/-2.6661e+00, -1.1187e-02,     9.0527e-05, -2.6430e-07,    3.2814e-10 /),    (/ 2.7708e+01,  1.3086e-01,    -1.5568e-03,  1.1635e-05,   -4.3171e-08,  6.8248e-11 /),   (/ 1.62180e-1, -1.16909e-3,     5.55185e-6, -8.75443e-9,     9.91567e-12,  0.00e+00   /) ] 
//       
         GasParV Water =    ["H2O   ", 1.310e1, 6.4714e+02,  2.2048320e+07,  5.6000e-02,  2.290e-01, 3.4400e-01,  1.8015000e-02,  1.8000e-00,  4.395e+00, -4.1860e-03,  1.4050000e-05, -1.5640e-08,  0.632e-11, 0.0000e+00,  0.0000000e+00,  0.00000e00,  0.0e0, (/ 0.00e+00,    0.00e+00,       0.00e+00,    0.000e+00, 0.00e+00 /), (/ 0.00e+00,    0.00e+00,       0.00e+00,    0.000e+00, 0.00e+00,    0.00e+00 /), (/ 0.00e+00,    0.00e+00,       0.00e+00,    0.000e+00, 0.00e+00,    0.00e+00 /) ]
//        
         GasParV Ethanol =  ["C2H5OH", 5.177e1, 5.1264e+02,  6.383475e+06,   1.670e-01,   2.480e-01, 5.5650e-01,  3.204200e-02,   1.700e-00,   4.714e+00, -6.9860e-03,  4.211000e-05,  -4.443e-08,   1.535e-11, 0.0000e+00,  0.0000000e+00,  0.00000e0,   0.0e0, (/ 0.00e+00,    0.00e+00,       0.00e+00,    0.000e+00,  0.00e+00 /), (/ 0.00e+00,    0.00e+00,       0.00e+00,    0.000e+00,   0.00e+00,    0.00e+00 /), (/ 0.00e+00,    0.00e+00,       0.00e+00,    0.000e+00, 0.00e+00,    0.00e+00 /) ] 
// 
// ----------
// ...... Derived-type parameter arrays     
// ----------
// 
*/
		 
		 
		 
		 //		 GasParV GasR[9] =  {  Methane,  Ethane,  Propane,  H2Sulfide,  CarbonO2, Nitrogen, Oxygen,  Water,    Ethanol };


		 //GasParV GasR[9];

		 GasParV GasR[9] =  {  Methane,  Ethane,  Propane,  H2Sulfide,  CarbonO2, Nitrogen, Oxygen,  Water,    Ethanol };


// 
// ----------
// ...... Derived-type arrays     
// ----------
// 
         GasPar2 Gas[9];
// 
// ----------
// ...... Double precision arrays
// ----------
// 

		 
		 /*
         float RK[NCom],Y[NCom],TR[NCom],Phi_F[NCom],ActCo[NCom];
         float AxAL[NCom],BB[NCom],ALF[NCom];
//
         float Kij_SC[NCom][NCom],Kij_PC[NCom][NCom],Kij_RC[NCom][NCom];
         float Kij_S[NCom][NCom],Kij_P[NCom][NCom],Kij_R[NCom][NCom];
         float G_1p[NCom][NCom],G_2p[NCom][NCom],Ex_a[NCom][NCom];

		 */





		          
		 float RK[9],Y[9],TR[9],Phi_F[9],ActCo[9];
         float AxAL[9],BB[9],ALF[9];
//
         float Kij_SC[9][9],Kij_PC[9][9],Kij_RC[9][9];
         float Kij_S[9][9],Kij_P[9][9],Kij_R[9][9];
         float G_1p[9][9],G_2p[9][9],Ex_a[9][9];
// 
// ----------
// ...... Integer arrays
// ----------
// 
         int ID[9];
// 
// ----------
// ...... Double precision variables
// ----------
// 
         float AUP,BUP,AxAL_m;
// 
// ----------
// ...... Integer variables
// ----------
// 
         int NumCom,iH2O,iH2S;
// 
// ----------
// ...... Character variables
// ----------
// 
         string EOS_Type;
         char Phase;
// 
// 
//      END MODULE RefRGas_Param



//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//








		 
//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//






		 
//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//




//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//

void Allocate_MemGen();
void PRINT_Headin();
void InOutFiles();
void Allocate_Input();
void READ_MainInputF();
void READ_RockData();
void READ_HydrateData();
void READ_ParamData();
void READ_ElemData();
void READ_ConxData();
void READ_GnerData();
void READ_InConData();
void NumFloPo_Digits();

void RGasSetup(int N_Com,string Cname[],string EquTyp,char PhaseF);
void SIMULATION_Cycle();


int main ()
{

//===================================== PROPERTIES OF AVAILABE SUBSTANCES ============================


//...... NOTE: The Psi, A, B in the parameters
// ......    For CO2 : From Diamond and Atkinfiev (2003), Fluid Phase Equilib (in press) 
// ......    For rest: Table 1 and 2, Atkinfiev and Diamond, 2003, GCA (see above)
// 
// ...... NOTE: Salting-out coefficients
// ......    Langmuir (1997)  Aqueous Environmental Geochemistry, Prentice Hall, Table 4.5, p.144 
// ......    Modified because of s problem with the log form in the reference



//============================ METHANE PROPERTIES

	Methane.Name=   "CH4   "     ;
	Methane.AtomDV= 2.514e1      ;
	Methane.TCrit=  1.9056e+02   ;
	Methane.PCrit=  4.600155e+06 ;
	Methane.VCrit=  9.9000e-02   ;
	Methane.ZCrit=  2.880e-01    ;
	Methane.Omega=  1.1000e-02   ;
	Methane.MolWt=  1.604300e-02 ;
	Methane.DMom=   0.0000e+00   ;
	Methane.A0=     4.568e+00    ;
    Methane.A1=     -8.9750e-03  ;
	Methane.A2=     3.631000e-05 ;
	Methane.A3=     -3.4070e-08  ;
	Methane.A4=     1.091e-11    ; 
	Methane.Psi=    -0.1131e+00  ;
	Methane.A=      -1.184620e+01;
	Methane.B=      1.48615e01   ;
	Methane.T_B=    111.6e0      ;
//	Methane.LogK=   [-2.5764e+00, -1.3019e-02,    1.0924e-04,-2.9008e-07,3.2774e-10];
	
	Methane.LogK[0]= -2.5764e+00 ;
	Methane.LogK[1]= -1.3019e-02 ;
	Methane.LogK[2]= 1.0924e-04  ;
	Methane.LogK[3]= -2.9008e-07 ;
	Methane.LogK[4]= 3.2774e-10  ;

	Methane.MVol[0]=3.7076e+01   ;
	Methane.MVol[1]=2.1313e-03   ;
	Methane.MVol[2]=-2.8752e-05  ;
	Methane.MVol[3]=1.1936e-06   ;
	Methane.MVol[4]=-7.3216e-09  ;
	Methane.MVol[5]=1.8923e-11   ;

	Methane.SOCo[0]=1.64818e-1  ;
	Methane.SOCo[1]=-1.40166e-3 ;
	Methane.SOCo[2]=1.32360e-5  ;
	Methane.SOCo[3]=-4.85733e-8 ;
	Methane.SOCo[4]=7.87967e-11 ;
	Methane.SOCo[5]=-5.52586e-14;



//    Alternative for CH4:  25-300 C 350 bar: Regression from Hnedkovsky, MVOL = ...
//    (/ 3.5821E+01, 5.8603E-02,-1.3091E-03, 1.9022E-05,-8.6735E-08, 1.4126E-10 /) 
//

	//============================ ETHANE PROPERTIES

	Ethane.Name=   "C2H6  "     ;
	Ethane.AtomDV= 4.566e1      ;
	Ethane.TCrit=  3.0532e+02   ;
	Ethane.PCrit=  4.883865e+06 ;
	Ethane.VCrit=  1.4800e-01   ;
	Ethane.ZCrit=  2.850e-01    ;
	Ethane.Omega=  9.9000e-02   ;
	Ethane.MolWt=  3.007000e-02 ;
	Ethane.DMom=   0.0000e+00   ;
	Ethane.A0=     4.178e+00    ;
    Ethane.A1=     -4.4270e-03  ;
	Ethane.A2=     5.660000e-05 ;
	Ethane.A3=     -6.6510e-08  ;
	Ethane.A4=     2.487e-11    ;
	Ethane.Psi=    -0.6091e+00  ;
	Ethane.A=      -1.634820e+01;
	Ethane.B=      2.00628e01   ;
	Ethane.T_B=    184.6e0      ;
	
	Ethane.LogK[0]= -2.1993e+00;
	Ethane.LogK[1]= -1.8015e-02;
	Ethane.LogK[2]=  1.4569e-04;
	Ethane.LogK[3]= -4.0442e-07;
	Ethane.LogK[4]=  4.6619e-10;

	Ethane.MVol[0]=   4.6452e+01 ;
	Ethane.MVol[1]=   2.3782e-01 ;
	Ethane.MVol[2]=   -2.8260e-03;
	Ethane.MVol[3]=   2.0155e-05 ;
	Ethane.MVol[4]=   -7.1926e-08;
	Ethane.MVol[5]=   1.0645e-10;


	Ethane.SOCo[0]= 1.64818e-1   ;
	Ethane.SOCo[1]=  -1.40166e-3 ;
	Ethane.SOCo[2]= 1.32360e-5   ;
	Ethane.SOCo[3]=  -4.85733e-8 ;
	Ethane.SOCo[4]= 7.87967e-11  ;
	Ethane.SOCo[5]=  -5.52586e-14;
	//============================ PROPANE PROPERTIES
  
	Propane.Name=   "C3H8  "      ;
	Propane.AtomDV= 6.618e1       ;
	Propane.TCrit=  3.6983e+02    ;
	Propane.PCrit=  4.2455175e+06 ;
	Propane.VCrit=  2.0300e-01    ;
	Propane.ZCrit=  2.810e-01     ;
	Propane.Omega=  1.5200e-01    ;  
	Propane.MolWt=  4.4097000e-02 ;
	Propane.DMom=   0.0000e+00    ;
	Propane.A0=     3.487e+00     ;  
    Propane.A1=     5.1310e-03    ;  
	Propane.A2=     6.0110000e-05 ; 
	Propane.A3=     -7.8930e-08   ; 
	Propane.A4=     3.079e-11     ;
	Propane.Psi=    -1.1471e+00   ;
	Propane.A=      -2.5387900e+01; 
	Propane.B=      2.82616e01    ;  
	Propane.T_B=    231.1e0       ;
	
	Propane.LogK[0]= -2.3150e+00;
	Propane.LogK[1]= -2.0996e-02;
	Propane.LogK[2]=  1.7215e-04;
	Propane.LogK[3]= -4.8268e-07;
	Propane.LogK[4]=  2.0620e-10;


	Propane.MVol[0]= 6.1525e+01;
	Propane.MVol[1]= 2.8840e-01;
	Propane.MVol[2]=-3.4270e-03;
	Propane.MVol[3]= 2.4452e-05;
	Propane.MVol[4]= -8.7292e-08;
	Propane.MVol[5]=1.2927e-10;

	Propane.SOCo[0]=0.00e+00;
	Propane.SOCo[1]=0.00e+00;
	Propane.SOCo[2]=0.00e+00;
	Propane.SOCo[3]=0.000e+00;
	Propane.SOCo[4]=0.00e+00;
	Propane.SOCo[5]=0.00e+00;


	//============================ H2S PROPERTIES		 


	H2Sulfide.Name="H2S   ";
	H2Sulfide.AtomDV=2.752e1;
	H2Sulfide.TCrit=3.7340e+02;
	H2Sulfide.PCrit=8.9368650e+06;
	H2Sulfide.VCrit=9.8500e-02;
	H2Sulfide.ZCrit=2.840e-01;
	H2Sulfide.Omega=9.0000e-02;
	H2Sulfide.MolWt=3.4082000e-02;
	H2Sulfide.DMom=9.0000e-01;
	H2Sulfide.A0=3.259e+00;
    H2Sulfide.A1=-3.4380e-03;
	H2Sulfide.A2=1.3190000e-05;
	H2Sulfide.A3=-1.3310e-08;
	H2Sulfide.A4=0.488e-11;
	H2Sulfide.Psi=-0.2102e+00;
	H2Sulfide.A= -1.1230300e+01;
	H2Sulfide.B=1.26104e01;
	H2Sulfide.T_B=0.0e0;
	
	H2Sulfide.LogK[0]=-6.7372e-01;
	H2Sulfide.LogK[1]=-1.4245e-02;
	H2Sulfide.LogK[2]=7.9816e-05;
	H2Sulfide.LogK[3]=-1.9861e-07;
	H2Sulfide.LogK[4]= 5.5931e-10;
	
	H2Sulfide.MVol[0]=3.3065e+01;
	H2Sulfide.MVol[1]=8.3733e-02;
	H2Sulfide.MVol[2]=-9.9447e-04;
	H2Sulfide.MVol[3]=6.9408e-06;
	H2Sulfide.MVol[4]=-2.4299e-08;
	H2Sulfide.MVol[5]=3.4721e-11;


	H2Sulfide.SOCo[0]=2.9050e-1;
	H2Sulfide.SOCo[1]=-1.57400e-4;
	H2Sulfide.SOCo[2]=-4.62000e+1;
	H2Sulfide.SOCo[3]=5.70500e-1;
	H2Sulfide.SOCo[4]=-1.7770e-3;
	H2Sulfide.SOCo[5]=0.00e+00;



// 
//    Alternative for H2S:  25-350 C 350 bar: Regression from Hnedkovsky, MVOL = ...
//    (/ 3.3521E+01, 6.8759E-02,-7.8118E-04, 9.2859E-06,-3.8445E-08, 6.1362E-11 /) 
// 

//============================ CO2 PROPERTIES
	
	CarbonO2.Name=   "CO2   "      ;
	CarbonO2.AtomDV= 2.690e1       ;
	CarbonO2.TCrit=  3.0412e+02    ;
	CarbonO2.PCrit=  7.3764600e+06 ;
	CarbonO2.VCrit=  9.4000e-02    ;
	CarbonO2.ZCrit=  2.740e-01     ;
	CarbonO2.Omega=  2.2500e-01    ;
	CarbonO2.MolWt=  4.4010000e-02 ;
	CarbonO2.DMom=   0.0000e-00    ;
	CarbonO2.A0=     3.539e+00     ;
    CarbonO2.A1=     1.3560e-03    ;
	CarbonO2.A2=     1.5020000e-05 ;
	CarbonO2.A3=     -2.3740e-08   ;
	CarbonO2.A4=     1.056e-11     ;
	CarbonO2.Psi=    -0.0880e+00   ;
	CarbonO2.A=      -9.3134000e+00;
	CarbonO2.B=      1.15477e01    ;
	CarbonO2.T_B=    0.0e0         ;
	
	CarbonO2.LogK[0]= -1.1161e+00;
	CarbonO2.LogK[1]= -1.6023e-02;
	CarbonO2.LogK[2]= 9.5730e-05;
	CarbonO2.LogK[3]= -2.3733e-07;
	CarbonO2.LogK[4]= 2.3441e-10;

	CarbonO2.MVol[0]= 2.9844e+01;
	CarbonO2.MVol[1]= 1.3722e-01;
	CarbonO2.MVol[2]=-1.6281e-03;
	CarbonO2.MVol[3]= 1.0898e-05;
	CarbonO2.MVol[4]= -3.6679e-08;
	CarbonO2.MVol[5]=  4.8444e-11;


	CarbonO2.SOCo[0]= -1.0312e00;
	CarbonO2.SOCo[1]= 1.28060e-3;
	CarbonO2.SOCo[2]= 2.55900e+2;
	CarbonO2.SOCo[3]= 4.44500e-1;
	CarbonO2.SOCo[4]= -1.6060e-3;
	CarbonO2.SOCo[5]= 0.00e+00;




// 
//    Alternative 1 for CO2:  25-350 C 350 bar: Regression from Hnedkovsky, MVOL = ...
//    (/ 3.3101E+01, 6.1728E-03, 3.2205E-06, 6.3807E-06,-3.6252E-08, 6.8151E-11 /)
// 
//    Alternative 2 for CO2:  25-350 C 350 bar: Regression of J.E. Garcia, LBNL, 2001, MVOL = ...
//    (/ 3.7510E+01,-9.5850E-02,+8.7400E-04,-5.0440E-07,-4.4956E-19, 6.3527E-22 /) 

// 
//============================ NITROGEN PROPERTIES
	
	
	Nitrogen.Name=    "N2    "        ;
	Nitrogen.AtomDV=  1.850e1         ;
	Nitrogen.TCrit=   1.2620e+02      ;
	Nitrogen.PCrit=   3.3943875e+06   ;
	Nitrogen.VCrit=   8.9500e-02      ;
	Nitrogen.ZCrit=   2.900e-01       ;
	Nitrogen.Omega=   3.7000e-02      ;
	Nitrogen.MolWt=   2.8014000e-02   ;
	Nitrogen.DMom=    0.0000e-00      ;
	Nitrogen.A0=      3.539e+00       ;
    Nitrogen.A1=      -0.2610e-03     ;
	Nitrogen.A2=      0.0070000e-05   ;
	Nitrogen.A3=      0.1570e-08      ;
	Nitrogen.A4=      -0.099e-11      ;
	Nitrogen.Psi=     -0.0320e+00     ;
	Nitrogen.A=       -1.1538000e+01  ;
	Nitrogen.B=       1.46278e01      ;
	Nitrogen.T_B=     0.0e0           ;
	
	Nitrogen.LogK[0]=-2.9831e+00;
	Nitrogen.LogK[1]=-9.9098e-03;
	Nitrogen.LogK[2]=8.5586e-05;
	Nitrogen.LogK[3]=-2.4931e-07;
	Nitrogen.LogK[4]=3.0552e-10;

	Nitrogen.MVol[0]= 3.0349e+01;
	Nitrogen.MVol[1]= 1.3953e-01;
	Nitrogen.MVol[2]=-1.6592e-03;
	Nitrogen.MVol[3]=1.2163e-05;
	Nitrogen.MVol[4]=-4.4426e-08;
	Nitrogen.MVol[5]= 6.8446e-11;

	Nitrogen.SOCo[0]= 8.68589e-2;
	Nitrogen.SOCo[1]=0.00e+00;
	Nitrogen.SOCo[2]=0.00e+00;
	Nitrogen.SOCo[3]=0.000e+00;
	Nitrogen.SOCo[4]=0.00e+00;
	Nitrogen.SOCo[5]=0.00e+00;



//============================ OXYGEN PROPERTIES

		
	Oxygen.Name=   "O2    "      ;
	Oxygen.AtomDV= 1.630e1       ;
	Oxygen.TCrit=  1.5458e+02    ;
	Oxygen.PCrit=  5.0459850e+06 ;
	Oxygen.VCrit=  7.3400e-02    ;
	Oxygen.ZCrit=  2.880e-01     ;
	Oxygen.Omega=  2.1000e-02    ;
	Oxygen.MolWt=  3.1999000e-02 ;
	Oxygen.DMom=   0.0000e-00    ;
	Oxygen.A0=     3.630e+00     ;
    Oxygen.A1=     -1.7940e-03   ;
	Oxygen.A2=     0.6580000e-05 ;
	Oxygen.A3=     -0.6010e-08   ;
	Oxygen.A4=     0.179e-11     ;
	Oxygen.Psi=    0.0260e+00    ;
	Oxygen.A=      -9.7540000e+00;
	Oxygen.B=      0.00000e00    ;
	Oxygen.T_B=    0.0e0         ;
	
	Oxygen.LogK[0]= -2.6661e+00;
	Oxygen.LogK[1]= -1.1187e-02;
	Oxygen.LogK[2]=  9.0527e-05;
	Oxygen.LogK[3]= -2.6430e-07;
	Oxygen.LogK[4]=  3.2814e-10;


	Oxygen.MVol[0]= 2.7708e+01 ;
	Oxygen.MVol[1]= 1.3086e-01 ;
	Oxygen.MVol[2]= -1.5568e-03;
	Oxygen.MVol[3]= 1.1635e-05 ;
	Oxygen.MVol[4]=-4.3171e-08 ;
    Oxygen.MVol[5]= 6.8248e-11 ;

	Oxygen.SOCo[0]= 1.62180e-1 ;
	Oxygen.SOCo[1]=-1.16909e-3 ;
	Oxygen.SOCo[2]= 5.55185e-6 ;
	Oxygen.SOCo[3]=-8.75443e-9 ;
	Oxygen.SOCo[4]= 9.91567e-12;
	Oxygen.SOCo[5]= 0.00e+00   ;



//============================ WATER PROPERTIES

	Water.Name=   "H2O   "     ;
	Water.AtomDV= 1.310e1      ;
	Water.TCrit=  6.4714e+02   ;
	Water.PCrit=  2.2048320e+07;
	Water.VCrit=  5.6000e-02   ;
	Water.ZCrit=  2.290e-01    ;
	Water.Omega=  3.4400e-01   ;
	Water.MolWt=  1.8015000e-02;
	Water.DMom=   1.8000e-00   ;
	Water.A0=     4.395e+00    ;
    Water.A1=     -4.1860e-03  ;
	Water.A2=     1.4050000e-05;
	Water.A3=     -1.5640e-08  ;
	Water.A4=     0.632e-11    ;
	Water.Psi=    0.0000e+00   ;
	Water.A=      0.0000000e+00;
	Water.B=      0.00000e00   ;
	Water.T_B=    0.0e0        ;
	
	Water.LogK[0]= 0.00e+00 ;
	Water.LogK[1]= 0.00e+00 ;
	Water.LogK[2]= 0.00e+00 ;
	Water.LogK[3]= 0.000e+00;
	Water.LogK[4]= 0.00e+00 ;


	Water.MVol[0]= 0.00e+00 ;
	Water.MVol[1]= 0.00e+00 ;
	Water.MVol[2]= 0.00e+00 ;
	Water.MVol[3]= 0.000e+00;
	Water.MVol[4]= 0.00e+00 ;
	Water.MVol[5]= 0.00e+00 ;


	Water.SOCo[0]= 0.00e+00;
	Water.SOCo[1]= 0.00e+00;
	Water.SOCo[2]= 0.00e+00;
	Water.SOCo[3]= 0.000e+00;
	Water.SOCo[4]= 0.00e+00;
	Water.SOCo[5]=0.00e+00;



//============================ ETHANOL PROPERTIES


	Ethanol.Name=   "C2H5OH";
	Ethanol.AtomDV= 5.177e1;
	Ethanol.TCrit=  5.1264e+02 ;
	Ethanol.PCrit=  6.383475e+06   ;
	Ethanol.VCrit=  1.670e-01;
	Ethanol.ZCrit=  2.480e-01 ;
	Ethanol.Omega=  5.5650e-01;
	Ethanol.MolWt=  3.204200e-02   ;
	Ethanol.DMom=   1.700e-00 ;
	Ethanol.A0=     4.714e+00;
    Ethanol.A1=     -6.9860e-03;
	Ethanol.A2=     4.211000e-05;
	Ethanol.A3=     -4.443e-08;
	Ethanol.A4=     1.535e-11;
	Ethanol.Psi=    0.0000e+00;
	Ethanol.A=      0.0000000e+00;
	Ethanol.B=      0.00000e0;
	Ethanol.T_B=    0.0e0;
	
	
	Ethanol.LogK[0]= 0.00e+00;
	Ethanol.LogK[1]= 0.00e+00;
	Ethanol.LogK[2]= 0.00e+00;
	Ethanol.LogK[3]= 0.000e+00;
	Ethanol.LogK[4]= 0.00e+00;


	Ethanol.MVol[0]= 0.00e+00;
	Ethanol.MVol[1]= 0.00e+00;
	Ethanol.MVol[2]= 0.00e+00;
	Ethanol.MVol[3]= 0.000e+00;
	Ethanol.MVol[4]= 0.00e+00;
	Ethanol.MVol[5]= 0.00e+00;
	
	Ethanol.SOCo[0]= 0.00e+00;
	Ethanol.SOCo[1]= 0.00e+00;
	Ethanol.SOCo[2]= 0.00e+00;
	Ethanol.SOCo[3]= 0.000e+00;
	Ethanol.SOCo[4]= 0.00e+00;
	Ethanol.SOCo[5]= 0.00e+00;


// 
		 
	GasR[0]=Methane;
	GasR[1]=Ethane;
	GasR[2]=Propane;
	GasR[3]=H2Sulfide;
	GasR[4]=CarbonO2;
	GasR[5]=Nitrogen;
	GasR[6]=Oxygen;
	GasR[7]=Water;
	GasR[8]=Ethanol;
	
//===================================== END OF PROPERTIES  ============================




//...... Modules to be used 
// 
//         USE Basic_Param
//         USE GenControl_Param


	char VV[5];

// -------
//... Double Precision Variables    
//-------
	
	float   ELT,ELT1,ELTC;

// -------
// ... Integer Variables    
// -------
// 
	
	int MC,i,ierG;

    
	string TITLE;
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of ISHR
//
//
// 
// -------
// ... Open file VERS to record subroutine use
// -------


	//INQUIRE(FILE='VERS',EXIST = EX_FileVERS)

	ofstream inVERS ("VERS");

	if (inVERS.is_open())
	{
		EX_FileVERS=true;
		cout << "able to open file";
		cout<< endl;
	}
	else
	{
		cout << "Unable to open file" << endl;
		cout<< endl;

		EX_FileVERS=false;
		inVERS.open("VERS");
	}




    
	cout<< "________________________________________________________________________________";
	cout<< "________________________________________________________________________________";
	cout<< "________________________________________________________________________________";
	cout<< " ******* ******      ******* *******        ******     ********    ***********"<<endl;
	cout<< " ******* ******      ******* *********    ***    ***   **********  ***********"<<endl;
	cout<< "   ***    *****      *****   ***    ***  ***      ***  ***    ***  *   ***   *"<<endl;
	cout<< "   ***    ******    ******   ***    *** ***        *** ***    ***      ***"<<endl;
	cout<< "   ***    *** ***  *** ***   *********  ***        *** ********        ***"<<endl;
	cout<< "   ***    ***  ******  ***   *******    ***        *** *******         ***"<<endl;
	cout<< "   ***    ***   ****   ***   ***         ***      ***  ***   ***       ***"<<endl;
	cout<< " ******* *****   **   ***** *****         ***    ***  *****   ***     *****"<<endl;
	cout<< " ******* *****        ***** *****           ******    *****    ***    *****"<<endl;
	cout<< "________________________________________________________________________________";
	cout<< "________________________________________________________________________________";
	cout<< "________________________________________________________________________________";
		
	cin>> TITLE;
	
	cout<<endl;

	cout<<"Problem Title:"<< TITLE<<endl;;



//
//***********************************************************************
//*                                                                     *
//*                Allocate memory for most arrays                      *
//*                                                                     *
//***********************************************************************
// 


		

	Allocate_MemGen();

// 
// 
//***********************************************************************
//*                                                                     *
//*          Print headings and important general information           *
//*                                                                     *
//***********************************************************************
// 
// 
		

	PRINT_Headin();

// 
//***********************************************************************
//*                                                                     *
//*               Check the directory for various files                 *
//*                                                                     *
//***********************************************************************
// 
// 
      


	InOutFiles();



//***********************************************************************
//*                                                                     *
//*                       Initialize the clock                          *
//*                                                                     *
//***********************************************************************
// 
// 

//	CPU_TimeF(TZERO);


//***********************************************************************
////*                                                                     *
//*             Allocate memory to temporary input arrays               *
//*                                                                     *
//***********************************************************************
 
 
      
	Allocate_Input();

/*
***********************************************************************
*                                                                     *
*             Read data from the main input data file                 *
*                                                                     *
***********************************************************************
*/

      //READ_MainInputF(MC)
	  
	READ_MainInputF();

	
	// 

// -------
//... Determine the floating point accuracy of the processor 
// -------



      
	
	NumFloPo_Digits();


//***********************************************************************
//*                                                                     *
//*      Read input data from the various files, both pre-existing      *
//*                  and created by the INPUT routine                   *
//*                                                                     *
//***********************************************************************
// 

	READ_Files();

// ***********************************************************************
// *                                                                     *
// *     Perform simulation by solving the mass and heat equations       *
// *                                                                     *
// ***********************************************************************
//


		SIMULATION_Cycle();


	  return 0;
}




//void EOS_DefaultNum (string EOS_Name, int MaxNum_MComp, int MaxNum_Equat, int MaxNum_Phase, int MaxNum_MobPh, int M_Add)
void EOS_DefaultNum (string EOS_Name)
{
//	int MaxNum_MComp,MaxNum_Equat,  MaxNum_Phase,MaxNum_MobPh,M_Add;

// -------------
// ...... Character variables
// -------------


//	char EOS_name[15];
//
//
// =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of EqTemp_HydrCH4 
//
//


	cout<< "EOS_DefaultNum    1.0   29 September 2004',6X, 'Default parameters of the hydrate equation of state"<<endl;

//
//
//***********************************************************************
//*                                                                     *
//*         DEFAULT NUMBERS FOR THE HYDRATE EQUATION OF STATE           *
//*                                                                     *
//***********************************************************************
// 
// 
    
	if(EOS_Name == "HYDRATE-EQU")
	{
		MaxNum_MComp = 3;   // Maximum number of mass components
		MaxNum_Equat = 4;   // Maximum number of equations
		MaxNum_Phase = 4;   // Maximum number of phases
		MaxNum_MobPh = 2;   // Maximum number of mobile phases
	}
	else if(EOS_Name == "HYDRATE-KIN")
	{
		MaxNum_MComp = 4;  // Maximum number of mass components
		MaxNum_Equat = 5;  // Maximum number of equations
		MaxNum_Phase = 4;  // Maximum number of phases
		MaxNum_MobPh = 2;  // Maximum number of mobile phases
	}
	else
	{
		cout<<"ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR -";
		cout<<endl;
		cout<< "                                 S I M U L A T I O N   A B O R T E D"<<endl;
		cout<< "The equation of state specified by EOS_Name  is unavailable"<<endl;
		cout<< "                                         PLEASE CORRECT AND TRY AGAIN";
		cin.get();
		cin.get();

		exit(0);
	}

	   //
// ----------
// ...... Adjust the parameter controlling the size of the auxilliary parameter array
// ----------
//

	if((EOS_Name == "HYDRATE") && (M_Add < 3)) 
	{
		M_Add = 3;
	}



}



void Allocate_MemGen()
{
//#include "EOS_Parameters.h"
//#include "EOS_DefaultParam.h"
//#include "Basic_Param.h"
//#include "GenControl_Param.h"
//#include "Diffusion_Param.h"
//#include "Element_Arrays.h"
//#include "Connection_Arrays.h"
//#include "Variable_Arrays.h"
//#include "SolMatrix_Arrays.h"
//#include "Q_Arrays.h"
//#include "PFMedProp_Arrays.h"
// 
//***********************************************************************
//***********************************************************************
//*                                                                     *
//*     MAIN ROUTINE FOR ALLOCATING MEMORY TO MOST ISHR ARRAYS          *
//*                                                                     *
//***********************************************************************
//***********************************************************************
	//IMPLICIT NONE
// 
// ----------
// ...... Integer Arrays    
// ----------
// 
	int ierror[24] = {999,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999,999};

// 
// ----------
// ...... Integer variables     
// ----------
//
	int NREDM_t,MPRIM_t,ii;
	int i,j, iiii;
// 
// ----------
// ...... Character variables     
// ----------
// 
	//char B_name = "               ";
	string B_name = "               ";
	//char HEADR[26];

	string HEADR;

	// 
// ----------
// ...... Logical variables     
// ----------
// 
	bool EX;
// 
//         SAVE //

//
// =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of Allocate_MemGen
	cout<<endl;
	cout<<"ASHR_Allocate_MemGen   1.0   1 APR 2014 Allocate memory to most arrays based on size provided by input parameters";
	cout<<endl;

	//
//

	ifstream inALLOC ("ALLOC");
	if (inALLOC.is_open())
	{

	}
    else
	{
		inALLOC.close();
		ofstream outALLOC ("ALLOC");
//		outALLOC.close ("ALLOC");
//		ifstream inALLOC ("ALLOC");
	}
  


/*         INQUIRE(FILE='ALLOC',EXIST=EX)
         IF(EX) THEN
            OPEN(50,FILE='ALLOC',STATUS='OLD')
            REWIND 50
         ELSE      
            OPEN(50,FILE='ALLOC',STATUS='NEW')
         END IF */    

// 
// 
//***********************************************************************
//*                                                                     *
//*              READ NUMBERS OF ELEMENTS AND CONNECTIONS               *
//*                                                                     *
//***********************************************************************
// 
//

	cin>> HEADR; // Read header from input
//
    if (HEADR == "ISHR_MEMORY_ALLOCATION")
	{
		cout<<endl<<endl<< "*********************      MEMORY ALLOCATION STARTED      **********************";
		cout<< endl;
	}
	else
	{
		cout<<"ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR -              '       S I M U L A T I O N   A B O R T E D      '                    The header 'ISHR MEMORY ALLOCATION' is missing at the top of the input file" << endl << "                       CORRECT AND TRY AGAIN                    ";
		cout<< endl;
		cin.get();
		cin.get();
		exit(0);
	}


//
  
	cin >> B_name;
    //EOS_Name = TRIM(ADJUSTL(B_name));
	EOS_Name = B_name;

//
//
//***********************************************************************
//*                                                                     *
//*  READ THE NUMBER OF COMPONENTS, EQUATIONS, PHASES, SEC. VARIABLES   *
//*                                                                     *
//***********************************************************************
		 

 // 
//

	cin >> NK>>NEQ>>NPH>>M_BinDif>>M_add;
//
// ...... Number of supplementary parameters 
//

	if(M_add <= 0)
	{
		M_Add = 1 ;           
	};




// 
// 
//***********************************************************************
//*                                                                     *
//*   DETERMINE DEFAULT PARAMETERS FOR THE VARIOUS EQUATIONS OF STATE   *
//*                                                                     *
//***********************************************************************
// 
// 

// 
	
	//EOS_DefaultNum(EOS_Name,  MaxNum_MComp,MaxNum_Equat,   MaxNum_Phase,MaxNum_MobPh,M_Add);
	EOS_DefaultNum(EOS_Name);

// ...... Determine the number of mobile phases 
//
         
	Num_MobPhs = MaxNum_MobPh;
	NPH        = MaxNum_Phase;

		 // ...... For hydrates, ensure that NEQ = NK+1 
//

	if(EOS_Name == "HYDRATE" && NEQ != (NK+1))
	{
		cout<< "//,20('ERROR-'),//,T33,  '       S I M U L A T I O N   A B O R T E D',//,   T08,'The number of equation NEQ =',i2,     ' <= the number of mass components NK =',i2,     ' : In hydrate simulation, NEQ = NK+1 '   /,T32,'               CORRECT AND TRY AGAIN',/,    //,20('ERROR-'))";
		cout<< NEQ,NK;
		exit(0);
	}
//
// ----------
// ...... Ensure that NK,NEQ,NPH do not exceed the maximum numbers 
// ----------
//

	if (NK  > MaxNum_MComp || NEQ > MaxNum_Equat ||  NPH > MaxNum_Phase)
	{
		cout<<"//,20('ERROR-'),//,T33,  '       S I M U L A T I O N   A B O R T E D',//,     T10,'One or more of the parameters NK,NEQ,NPH = (',  i2,',',i2,',',i2,')',/,     T10,'exceed the maximum values ',      '"<<MaxNum_MComp<<", "<<MaxNum_Equat<<", ',       '"<<MaxNum_Phase<<" = (',i2,',',i2,',',i2,')',/,  /,T32,'               CORRECT AND TRY AGAIN',        //,20('ERROR-'))";
		cout<< NK,NEQ,NPH,    MaxNum_MComp,MaxNum_Equat,MaxNum_Phase ;
		exit(0);
	}

//
// ...... Number of secondary parameters (to be stored in the DP array)
//
		 
	if(M_BinDif == 0) 
	{
		Num_SecPar = 8*MaxNum_Phase+1+M_add;
	}
	else
	{
		Num_SecPar = 10*MaxNum_Phase+1+M_add;
	}
//
// ----------
// ...... Read the number of elements, connections and 
// ...... number of characters in the element name
// ----------
//

	cin>> MaxNum_Elem>>MaxNum_Conx>>No_CEN>>FLG_con;
//
// ...... NCONUP: The total number of connections of elements with 
// ......         other elements of larger element numbers 
// ...... NCONDN: The total number of connections of elements with
// ......         other elements of smaller element numbers 
//

	NCONUP = 11*MaxNum_Conx/10;
	NCONDN = NCONUP;

//
//
//***********************************************************************
//*                                                                     *
//*            MEMORY ALLOCATION TO ELEMENT-RELATED ARRAYS              *
//*                                                                     *
//***********************************************************************
//
//

	if(No_CEN != 8)
	{
		No_CEN = 5;
	}

// 
// ----------
// ...... Allocate memory to derived-type arrays  
// ----------
// 
	
	elem= new Element [MaxNum_Elem];    // ALLOCATE(elem(MaxNum_Elem), STAT=ierror(1))


	
	if (!elem)
	{
		ierror[0]=1;
	}
	else
	{
		ierror[0]=0;
	}


// 
// ----------
// ...... Allocate memory to double precision allocatable arrays     
// ----------
// 

	AI= new float[MaxNum_Elem];    // ALLOCATE(AI(MaxNum_Elem),      STAT=ierror(2))
	if (!AI)
	{
		ierror[1]=1;
	}
	else
	{
		ierror[1]=0;
	}


	AHT= new float[MaxNum_Elem];    // ALLOCATE(AHT(MaxNum_Elem),     STAT=ierror(3))
	if (!AHT)
	{
		ierror[2]=1;
	}
	else
	{
		ierror[2]=0;
	}
// 

	Pm= new float[MaxNum_Elem];    // ALLOCATE(Pm(MaxNum_Elem),      STAT=ierror(4))
	if (!Pm)
	{
		ierror[3]=1;
	}
	else
	{
		ierror[3]=0;
	}


	X_Coord= new float[MaxNum_Elem];    // ALLOCATE(X_Coord(MaxNum_Elem),    Y_Coord(MaxNum_Elem),    Z_Coord(MaxNum_Elem), STAT=ierror(5))
	Y_Coord= new float[MaxNum_Elem];
	Z_Coord= new float[MaxNum_Elem];

	
		 
	if (!X_Coord || !Y_Coord || !Z_Coord)
	{
		ierror[4]=1;
	}
	else
	{
		ierror[4]=0;
	}
//
//
//***********************************************************************
//*                                                                     *
//*          MEMORY ALLOCATION TO CONNECTION-RELATED ARRAYS             *
//*                                                                     *
//***********************************************************************
// 
// ----------

// 
// ----------
// ...... Allocate memory to derived-type arrays  
// ----------
// 

	conx= new Connection[MaxNum_Conx]; //ALLOCATE(conx(MaxNum_Conx), STAT=ierror(6))
	//cout<< conx[0].area << conx[29].area<<conx[50].area; 
	if (!conx)
	{
		ierror[5]=1;
	}
	else
	{
		ierror[5]=0;
	}

// 

	if(MaxNum_Phase > 1) 
	{
		conV= new VaporFlux[MaxNum_Conx]; //ALLOCATE(conV(MaxNum_Conx), STAT=ierror(7))
	}
	if (!conV)
	{
		ierror[6]=1;
	}
	else
	{
		ierror[6]=0;
	}
		 // 
// ----------
// ...... Allocate memory to arrays within the derived type: Conx  
// ----------
// 

	
	/*

	for (i = 1; i<=MaxNum_Conx;1)
	{
			 //ALLOCATE(Conx(i).FluxF(Num_MobPhs),    Conx(i).VelPh(Num_MobPhs),     STAT = ierror(8))
		conx(i)->FluxF= new float[Num_MobPhs];
		conx(i)->VelPh= new float[Num_MobPhs];
	}
*/
	
	//conx->FluxF = new float[Num_MobPhs];
	//conx->VelPh = new float[Num_MobPhs];

	for (i=0;i<MaxNum_Conx;i++)
	{
		conx[i].FluxF = new float[Num_MobPhs];
		conx[i].VelPh = new float[Num_MobPhs];
	}

	if (!conx->FluxF || !conx->VelPh)
	{
		ierror[7]=1;
	}
	else
	{
		ierror[7]=0;
	}
// 
// ----------
// ...... Allocate memory to double precision allocatable arrays     
// ----------
// 

	sig= new float[MaxNum_Conx];   //ALLOCATE(sig(MaxNum_Conx), STAT=ierror(9))
	if (!sig) 
	{
		ierror[8]=1;
	}
	else
	{
		ierror[8]=0;
	}
//
//

//
//
//***********************************************************************
//*                                                                     *
//*      MEMORY ALLOCATION TO PRIMARY & SECONDARY VARIABLE ARRAYS       *
//*                                                                     *
//***********************************************************************
//
//

	MaxNum_PVar = (MaxNum_MComp+1)*MaxNum_Elem;
	NREDM       =  NEQ*MaxNum_Elem;
// 
	if(FLG_con == "ACTIVE") 
	{
		NREDM_t = 1;
		MPRIM_t = 1;
	}
	else
	{
		NREDM_t = NEQ*MaxNum_Elem;
		MPRIM_t = MaxNum_PVar;
	}

// 
// ----------
// ...... Allocate memory to double precision allocatable arrays     
// ----------
// 

	X= new float [MaxNum_PVar];        //ALLOCATE(X(MaxNum_PVar),       STAT=ierror(10))
	if (!X)
	{
		ierror[9]=1;
	}
	else
	{
		ierror[9]=0;
	}


	DX= new float [MPRIM_t];           //ALLOCATE(DX(MPRIM_t),          STAT=ierror(11))
	if (!DX)
	{
		ierror[10]=1;
		
	}
	else
	{
		ierror[10]=0;
	}


	DELX= new float [MPRIM_t];         //ALLOCATE(DELX(MPRIM_t),        STAT=ierror(12))
	if (!DELX)
	{
		ierror[11]=1;
	}
	else
	{
		ierror[11]=0;
	}
//

	R= new float [NEQ*MaxNum_Elem+1];  //ALLOCATE(R(NEQ*MaxNum_Elem+1), STAT=ierror(13))
	if (!R)
	{
		ierror[12]=1;
	}
	else
	{
		ierror[12]=0;
	}
// 

	DOLD= new float [NREDM_t];         //ALLOCATE(DOLD(NREDM_t),        STAT=ierror(14))
	if (!DOLD)
	{
		ierror[13]=1;
	}
	else
	{
		ierror[13]=0;
	}
// 
// ----------

// 
// ----------
// ...... Allocate memory to derived-type arrays - Secondary variables  
// ----------
// 

	Cell_V=new SecParam[MaxNum_Elem,NEQ+1];//ALLOCATE(Cell_V(MaxNum_Elem,0:NEQ), STAT=ierror(15))
	if (!Cell_V)
	{
		ierror[14]=1;
	}
	else
	{
		ierror[14]=0;
	}
// 
// ----------
// ...... Allocate memory for arrays within the derived type  
// ----------
// 

	Cell_V->p_Satr = new float[MaxNum_Phase];
	Cell_V->p_KRel = new float[MaxNum_Phase];
	Cell_V->p_Visc = new float[MaxNum_Phase];
	Cell_V->p_Dens = new float[MaxNum_Phase];
	Cell_V->p_Enth = new float[MaxNum_Phase];
	Cell_V->p_PCap = new float [MaxNum_Phase];
	Cell_V->p_MasF = new float [MaxNum_MComp,MaxNum_Phase];
			
	if (!Cell_V->p_Satr || !Cell_V->p_KRel || ! Cell_V->p_Visc || !Cell_V->p_Dens || !Cell_V->p_Enth || ! Cell_V->p_PCap || !Cell_V->p_MasF)
	{
		ierror[15]=1;
	}
	else
	{
		ierror[15]=0;
	}
	/*

	for (j = 0; j<=NEQ;1)
	{
		for (i=1; i<= MaxNum_Elem;1)
		{
			Cell_V(i,j).p_Satr= new float[MaxNum_Phase];
			Cell_V(i,j).p_KRel = new float [MaxNum_Phase];
			Cell_V(i,j).p_Visc = new float [MaxNum_Phase];
			Cell_V(i,j).p_Dens = new float [MaxNum_Phase];
			Cell_V(i,j).p_Enth = new float [MaxNum_Phase];
			Cell_V(i,j).p_PCap = new float [MaxNum_Phase];
			Cell_V(i,j).p_MasF = new float [MaxNum_MComp,MaxNum_Phase];
				 
			if (!Cell_V(i,j).p_Satr || !Cell_V(i,j).p_KRel || ! Cell_V(i,j).p_Visc || !Cell_V(i,j).p_Dens || !Cell_V(i,j).p_Enth || ! Cell_V(i,j).p_PCap || !Cell_V(i,j).p_MasF)
			{
				ierror[16]=1;
			}
			else
			{
				ierror[16]=0;
			}
		}
	}
*/




//         DO_j1 : DO j = 0,NEQ
 //           DO_i1 : DO i = 1,MaxNum_Elem
// 
 //              ALLOCATE(Cell_V(i,j)%p_Satr(MaxNum_Phase),  Cell_V(i,j)%p_KRel(MaxNum_Phase),     Cell_V(i,j)%p_Visc(MaxNum_Phase),     Cell_V(i,j)%p_Dens(MaxNum_Phase),     Cell_V(i,j)%p_Enth(MaxNum_Phase),       Cell_V(i,j)%p_PCap(MaxNum_Phase),       Cell_V(i,j)%p_MasF(MaxNum_MComp,MaxNum_Phase),   STAT = ierror(16))
// 
 //           END DO DO_i1
 //        END DO DO_j1
// 
// ----------
// ...... Allocate memory for diffusion coefficients in the derived-type arrays (M_BinDif /= 0)  
// ----------
// 
/*
	if(M_BinDif != 0 && MaxNum_MComp >= 2)
	{
		for (j = 0; j<=NEQ;1)
		{
			for (i = 0; i<=MaxNum_Elem;1)
			{
				Cell_V(i,j).p_DifA = new float[MaxNum_Phase];
				Cell_V(i,j).p_DifB = new float[MaxNum_Phase];  //ALLOCATE(Cell_V(i,j)%p_DifA(MaxNum_Phase),     Cell_V(i,j)%p_DifB(MaxNum_Phase),    STAT = ierror(17));

				if (!Cell_V(i,j).p_DifA || !Cell_V(i,j).p_DifB)
				{
					ierror[17]=1;
				}
				else
				{
					ierror[17]=0;
				}

			}
		} // end if


//            DO_j2 : DO j = 0,NEQ
//               DO_i2 : DO i = 1,MaxNum_Elem
//                  ALLOCATE(Cell_V(i,j)%p_DifA(MaxNum_Phase),     Cell_V(i,j)%p_DifB(MaxNum_Phase),    STAT = ierror(17))
//               END DO DO_i2
//            END DO DO_j2
// 

	}
*/

		
	if(M_BinDif != 0 && MaxNum_MComp >= 2)
	{
		Cell_V->p_DiFA= new float[MaxNum_Phase];
		Cell_V->p_DiFB = new float[MaxNum_Phase];  //ALLOCATE(Cell_V(i,j)%p_DifA(MaxNum_Phase),     Cell_V(i,j)%p_DifB(MaxNum_Phase),    STAT = ierror(17));
		if (!Cell_V->p_DiFA || !Cell_V->p_DiFB)
		{
			ierror[16]=1;
		}
		else
		{
			ierror[16]=0;
		}
	}





// 
// ----------
// ...... Allocate memory for additional data in the derived-type arrays 
// ----------
// 



//		 if(M_add > 0)
//		 {
//
//			 DO_j3 : DO j = 0,NEQ
//               DO_i3 : DO i = 1,MaxNum_Elem
//                 ALLOCATE(Cell_V(i,j)%p_AddD(M_add),   STAT = ierror(18))
//               END DO DO_i3
//            END DO DO_j3
// 
//         }
/*
	if(M_add > 0)
	{
		for (j = 0; j<=NEQ;1)
		{
			for (i=1; i<=MaxNum_Elem;1)
			{
				Cell_V(i,j).p_AddD = new float [M_add];

				if (//Cell_V(i,j).p_AddD)
				{
					ierror[18]=1;
				}
				else
				{
					ierror[18]=0;
				}
			}
		}
	}
*/

		if(M_add > 0)
		{
			Cell_V->p_AddD = new float [M_add];
			if (!Cell_V->p_AddD)
			{
				ierror[17]=1;
			}
			else
			{
				ierror[17]=0;
			}
		}


//
//
//***********************************************************************
//*                                                                     *
//*                  MEMORY ALLOCATION TO WORK ARRAYS                   *
//*                                                                     *
//***********************************************************************
//
//

	WKAREA= new float [NREDM+10];//ALLOCATE(WKAREA(NREDM+10), STAT=ierror(19))
	if (!WKAREA)
	{
		ierror[18]=1;
	}
	else
	{
		ierror[18]=0;
	}
//
//
//***********************************************************************
//*                                                                     *
//*          MEMORY ALLOCATION TO SOURCE/SINK-RELATED ARRAYS            *
//*                                                                     *
//***********************************************************************
//
//

	cin>>MaxNum_SS;

// 
// ----------
// ...... Allocate memory to derived-type     
// ----------

// 

	SS= new SourceSink[MaxNum_SS]; //ALLOCATE(SS(MaxNum_SS), STAT = ierror(20))

	if (!SS)
	{
		ierror[19]=1;
	}
	else
	{
		ierror[19]=0;
	}
// 

/*
	if(MaxNum_Phase > 1) 
	{
		for (i = 1; i<=MaxNum_SS;1)
		{
				 //ALLOCATE(SS(i)%frac_flo(Num_MobPhs),    SS(i)%rate_phs(Num_MobPhs),  STAT = ierror(21));

			SS(i).frac_flo = new float [Num_MobPhs];
			SS(i).rate_phs = new float [Num_MobPhs];

			if (//SS(i).frac_flo || //SS(i).rate_phs)
			{
				ierror[21]=1;
			}
			else
			{
				ierror[21]=0;
			}
		}
	}

*/



	if(MaxNum_Phase > 1) 
	{
		SS->frac_flo = new float [Num_MobPhs];
		SS->rate_phs = new float [Num_MobPhs];
		if (!SS->frac_flo || !SS->rate_phs)
		{
			ierror[20]=1;
		}
		else
		{
			ierror[20]=0;
		}
	}



	Well_TData= new Tabular[MaxNum_SS];//ALLOCATE(Well_TData(MaxNum_SS), STAT = ierror(21))
	if (!Well_TData)
	{
		ierror[20]=1;
	}
	else
	{
		ierror[20]=0;
	}
// 

	WDel = new Tab_BHP [MaxNum_SS];//ALLOCATE(WDel(MaxNum_SS),       STAT = ierror(22))
	if (!WDel)
	{
		ierror[21]=1;
	}
	else
	{
		ierror[21]=0;
	}

//
//
//***********************************************************************
//*                                                                     *
//*           MEMORY ALLOCATION TO ARRAYS OF ROCK PROPERTIES            *
//*                                                                     *
//***********************************************************************
//
//

	cin>> MaxNum_Media;
// 
// ----------
// ...... Allocate memory to derived-type allocatable arrays     
// ----------
// 

	PoMed= new PFMedium1 [MaxNum_Media];//ALLOCATE(PoMed(MaxNum_Media), STAT=ierror(23))

	if (!PoMed)
	{
		ierror[22]=1;
	}
	else
	{
		ierror[22]=0;
	}
// 

	media= new PFMedium [MaxNum_Media];//ALLOCATE(media(MaxNum_Media), STAT=ierror(24))

	if (!media)
	{
		ierror[23]=1;
	}
	else
	{
		ierror[23]=0;
	}


// 
// 
	
	
	for (iiii=0 ; iiii<=23 ; iiii++)
	{
		if(ierror[iiii] == 0)
		{
			cout<<"Memory allocation at point in subroutine Allocate_MemGen was successful  " << iiii;
			cout<< endl;

		}
		else if(ierror[iiii] == 999)
		{
			continue;
		}
		else
		{
			cout<< "ERROR- Memory allocation at point in subroutine Allocate_MemGen was unsuccessful//////////////";
			cout<<"THE EXECUTION WAS STOPPED    //////////////ERROR"<<iiii;
			cout<< endl;
			cin.get();
			cin.get();
			exit(0);
		}

	}

	cout<< "ALLOCATION FINISHED";


}




      
void PRINT_Headin()
{

// 
// ...... Modules to be used 
// 
         



//	#include "EOS_Parameters.h"
//	#include "Basic_Param.h"
 //   #include "GenControl_Param.h"


//C
//C***********************************************************************
//C***********************************************************************
//C*                                                                     *
//C*//         ROUTINE FOR PRINTING HEADINGS IN THE OUTPUT FILE           *
//C*//                   Version 1.0 - April 23, 2003                     *     
//C*                                                                     *
//C***********************************************************************
//C***********************************************************************
//C
//      IMPLICIT NONE
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of PRINT_Headin
//
//

//
// ... Format
//
	cout<< endl;
	cout<<"PRINT_Headin      1.0   1 April    2014";
	cout<< endl;
	cout<<"Routine for printing the headings in the output file";
	cout<< endl;


// 
// -------
// ... Print headers         
// -------
// 
    

	//cout<<"@   @  @   @  @@@   @@@    @@   @@@@@  @@@@     ',  '@@@   @@@@   @@@      @@@  @  @   @     @@@   @  @  @   @',/,  T14,'@   @  @   @  @  @  @  @  @  @    @    @        ', '@  @  @     @        @     @  @@ @@     @  @  @  @  @@  @',/,  T14,'@@@@@   @ @   @  @  @@@   @@@@    @    @@@@     ', '@@@   @@@@   @@       @@   @  @ @ @     @@@   @  @  @ @ @',/,  T14,'@   @    @    @  @  @ @   @  @    @    @        ', '@ @   @        @        @  @  @   @     @ @   @  @  @  @@',/,  T14,'@   @    @    @@@   @  @  @  @    @    @@@@.....', '@  @  @@@@  @@@ .....@@@   @  @   @     @  @   @@   @   @',//";
	cout<< endl;



//


	cout<<"ISHR IS A PROGRAM FOR MULTIPHASE MULTICOMPONENT FLOW IN PERMEABLE MEDIA, INCLUDING HEAT FLOW.";
	cout<<endl;
	cout<<"*****************************************************************************************************************************************************************************************************************************************";
	cout<< endl;

//
//
//

	cout<<"ISHR v1.0 IS WRITTEN IN C++. USING AN OBJECT_ORIENTED PROGRAMMING STRUCTURE AND WAS DEVELOPED BY B. GHAREDAGHLOO, MSC OF RESERVOIR ENGINEERING, KPE Co. IRAN";
	cout<< endl;




//
//
//
//
 
	cout<< "ISHR V1.0 (May 2014)";
	cout<< endl;


//

	//cout<< "HydrateResSim Copyright (c) 2004, The Regents', ' of the University of California (through Lawrence', ' Berkeley National Laboratory).',/, ' All rights reserved.')";

//
//
//

 
//	cout<<"//, ' NOTICE: Redistribution and use in source and binary forms,', ' with or without modification, are permitted provided that ',/, ' the following conditions are met:',//, ' (1) Redistributions of source code must retain the above',  ' copyright notice, this list of conditions and the ',/, '     following disclaimer.',/,  ' (2) Redistributions in binary form must reproduce the above', ' copyright notice, this list of conditions and the following',/,  '     disclaimer in the documentation and/or other materials', ' provided with the distribution.',/, ' (3) Neither the name of HydrateResSim, the National Energy', ' Technology Laboratory, the U.S. Dept. of Energy nor the', ' names',/, '     of its contributors may be used to endorse or promote', ' products derived from this software without specific',/, '     prior written permission.',//, ' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND',  ' CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED', ' WARRANTIES, INCLUDING',/,  ' BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY', ' AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.', ' IN NO EVENT',/, ' SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR', ' ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,', ' OR CONSEQUENTIAL',/, ' DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF', ' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;', ' OR BUSINESS',/, ' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,', ' WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING',/, ' NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE', ' OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY', ' OF SUCH DAMAGE.',/)";



//
      
	
	cout<< MaxNum_Elem<<" "<<MaxNum_Conx<<" "<<MaxNum_Equat<<" "<<MaxNum_MComp<<" "<<  MaxNum_Phase<<" "<<MaxNum_SS<<" "<<M_BinDif<<" "<<M_Add<<endl;
	//,' PARAMETERS FOR DYNAMIC DIMENSIONING OF MAJOR',   ' ARRAYS (MAIN PROGRAM) ARE AS FOLLOWS',//,   T10,'MaxNum_Elem  =',I8,/,T10,'MaxNum_Conx  =',I8,/,    T10,'MaxNum_Equat =',I8,/,T10,'MaxNum_MComP =',I8,/,   T10,'MaxNum_Phase =',I8,/,T10,'MaxNum_SS    =',I8,/,    T10,'M_BinDif     =',I8,/,T10,'M_Add        =',I8,//,  ' ',131('='))
//
//
//
      
	cout<< MaxNum_Elem<<" "<<MaxNum_Conx<<" "<<MaxNum_PVar<<" "<<   MaxNum_SS<<" "<<MaxNum_Media<<" "<<endl;

//
// ... Format
//
 
	  
	cout<<"MAXIMUM NUMBER OF VOLUME ELEMENTS (GRID BLOCKS):    "; 
	cout<<MaxNum_Elem<< endl;
	cout<<"MAXIMUM NUMBER OF CONNECTIONS (INTERFACES):         ";
	cout<<MaxNum_Conx<<endl;
	cout<<"MAXIMUM LENGTH OF PRIMARY VARIABLE ARRAYS:          ";
	cout<< MaxNum_PVar<<endl;
	cout<<"MAXIMUM NUMBER OF GENERATION ITEMS (SINKS/SOURCES): ";
	cout<<MaxNum_SS<<endl;
	cout<<"MAXIMUM NUMBER OF POROUS/FRACTURED MEDIA            ";
	cout<<MaxNum_Media<<endl;

//
//
//
      
	cout<<"Array dimensioning is dynamic. The use of the direct solver will reduce the size of the maximum solvable problem";
	cout<<endl;


//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End of PRINT_Headin
//

//
	  
}




void InOutFiles()
{

//         USE Basic_Param, ONLY: EX_FileVERS
//
// ***********************************************************************
// ***********************************************************************
//                                                                       *
// //  ROUTINE FOR CHECKING/OPENING THE FILES NEEDED BY "HydrateResSim"  *
//                                                                       *
////                  Version 1.0 - January 14, 2003                     *     
//                                                                       *
// ***********************************************************************
// ***********************************************************************
//
//      IMPLICIT NONE
// 
// ----------
// ... Logical variables
// ----------
// 
      bool EX;
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of InOutFiles
//
//
	  cout<<endl<<endl<<endl<<endl;
      cout<< "SUMMARY OF DISK FILES" <<endl;    // Write comment in the output file     
      cout<<"InOutFiles  1.0   14 May 2014, Open files VERS, MESH, INCON, GENER, SAVE LINEQ, and TABLE"<<endl;




      if(EX_FileVERS)
	  {
		  cout<<"FILE 'VERS'  EXISTS --- OPEN AS AN OLD FILE"<<endl;
	  }
	  else
	  {
		  
		  cout<<"FILE 'VERS'  DOES NOT EXIST --- OPEN AS A NEW FILE"<<endl;
		  
	  }

			   // 
// 
// -------
// ... File MESH
// -------
// 
      	

	  
	  
	  ifstream MESH_OUT ("MESH");
	  if (MESH_OUT.is_open())
	  {
		   cout<<"FILE 'MESH'  EXISTS --- OPEN AS AN OLD FILE"<<endl;
		   MESH_OUT.close();
		   MESH_OUT.open("MESH");
	  }
	  else
	  {
		  cout<<"FILE 'MESH'  DOES NOT EXIST --- OPEN AS A NEW FILE"<<endl;
		  MESH_OUT.close();
		  ofstream MESH_OUT ("MESH");
	  }
// 
// -------
// ... File INCON
// -------
// 
      ifstream INCON_OUT ("INCON");
	  if(INCON_OUT.is_open())
	  {
		   cout<<"FILE 'INCON'  EXISTS --- OPEN AS AN OLD FILE"<<endl;
		   INCON_OUT.close();
		   INCON_OUT.open("INCON");
	  }
	  else
	  {
		  cout<<"FILE 'INCON'  DOES NOT EXIST --- OPEN AS A NEW FILE"<<endl;
		  INCON_OUT.close();
		  ofstream INCON_OUT ("INCON");
	  }

// 
// -------
// ... File GENER
// -------
// 
      
	  ifstream GENER_OUT ("GENER");
	  if(GENER_OUT.is_open())
	  {
		   cout<<"FILE 'GENER'  EXISTS --- OPEN AS AN OLD FILE"<<endl;
		   GENER_OUT.close();
		   GENER_OUT.open("GENER");
	  }
	  else
	  {
		  cout<<"FILE 'GENER'  DOES NOT EXIST --- OPEN AS A NEW FILE"<<endl;
		  GENER_OUT.close();
		  ofstream GENER_OUT ("GENER");
	  }
	  
	  
//	  INQUIRE(FILE='GENER',EXIST = EX)
// 
//      IF_GENER: IF(EX) THEN
//                   PRINT 6008
//                   OPEN(3,FILE='GENER',STATUS='OLD')
//                ELSE
 //                  PRINT 6009
  //                 OPEN(3,FILE='GENER',STATUS='NEW')
  //                 ENDFILE 3
  //              END IF IF_GENER
// 
// -------
// ... File SAVE
// -------
//
	  ifstream SAVE_OUT ("SAVE");
	  if(SAVE_OUT.is_open())
	  {
		   cout<<"FILE 'SAVE'   EXISTS --- OPEN AS AN OLD FILE"<<endl;
		   SAVE_OUT.close();
		   SAVE_OUT.open("SAVE");
	  }
	  else
	  {
		  cout<<"FILE 'SAVE'   DOES NOT EXIST --- OPEN AS A NEW FILE"<<endl;
		  SAVE_OUT.close();
		  ofstream SAVE_OUT ("SAVE");
	  }


	  //      INQUIRE(FILE='SAVE',EXIST = EX)
// 
//      IF_SAVE: IF(EX) THEN
//                  PRINT 6010
//                  OPEN(7,FILE='SAVE',STATUS='OLD')
//               ELSE
//                  PRINT 6011
//                  OPEN(7,FILE='SAVE',STATUS='NEW')
//              END IF IF_SAVE
// 
// -------

	  
	  // ... File LINEQ
// -------
// 
	  ifstream LINEQ_OUT ("LINEQ");
	  if(LINEQ_OUT.is_open())
	  {
		   cout<<"FILE 'LINEQ'   EXISTS --- OPEN AS AN OLD FILE"<<endl;
		   LINEQ_OUT.close();
		   LINEQ_OUT.open("LINEQ");
	  }
	  else
	  {
		  cout<<"FILE 'LINEQ'   DOES NOT EXIST --- OPEN AS A NEW FILE"<<endl;
		  LINEQ_OUT.close();
		  ofstream LINEQ_OUT ("LINEQ");
	  }

 //     INQUIRE(FILE='LINEQ',EXIST = EX)
// 
//      IF_LINEQ: IF(EX) THEN
//                   PRINT 6012
//                   OPEN(15,FILE='LINEQ',STATUS='OLD')
//                ELSE
//                   PRINT 6013
//                   OPEN(15,FILE='LINEQ',STATUS='NEW')
//                   REWIND 15
//                END IF IF_LINEQ
// 
// -------
// ... File TABLE
// -------
// 

	  ifstream TABLE_OUT ("TABLE");
	  if(TABLE_OUT.is_open())
	  {
		   cout<<"FILE 'TABLE'   EXISTS --- OPEN AS AN OLD FILE"<<endl;
		   TABLE_OUT.close();
		   TABLE_OUT.open("TABLE");
	  }
	  else
	  {
		  cout<<"FILE 'TABLE'   DOES NOT EXIST --- OPEN AS A NEW FILE"<<endl;
		  TABLE_OUT.close();
		  ofstream TABLE_OUT ("TABLE");
	  }

	  
	  /*
	  INQUIRE(FILE='TABLE',EXIST = EX)
      IF_TABLE: IF(EX) THEN
                   PRINT 6014
                   OPEN(8,FILE='TABLE',STATUS='OLD')
                ELSE
                   PRINT 6015
                   OPEN(8,FILE='TABLE',STATUS='NEW')
                   ENDFILE 8
                END IF IF_TABLE
	*/


  
	  ifstream OUTPUTRESULT ("OUTPUTRESULT.out");
	  if(OUTPUTRESULT.is_open())
	  {
		   cout<<"FILE 'OUTPUTRESULT.out'   EXISTS --- OPEN AS AN OLD FILE"<<endl;
		   OUTPUTRESULT.close();
		   OUTPUTRESULT.open("OUTPUTRESULT.out");
	  }
	  else
	  {
		  cout<<"FILE 'OUTPUTRESULT.out'   DOES NOT EXIST --- OPEN AS A NEW FILE"<<endl;
		  OUTPUTRESULT.close();
		  ofstream OUTPUTRESULT ("OUTPUTRESULT.out");
	  }


      /*
	  INQUIRE(FILE='OUTPUTRESULT.out',EXIST = EX)

      IF_OUTPUTRESULT: IF(EX) THEN
                   PRINT 6016
                   OPEN(14,FILE='OUTPUTRESULT.out',STATUS='OLD')
                ELSE
                   PRINT 6017
                   OPEN(14,FILE='OUTPUTRESULT.out',STATUS='NEW')
                   ENDFILE 14
                
				END IF IF_OUTPUTRESULT
	
	
	   */


	  
}




/*
void CPU_TimeF(time)
{
	//IMPLICIT NONE


	//REAL(KIND = 8), INTENT(OUT) :: time

	float time;
	int itime,irate;

	int icall = 0;
//
   //   SAVE icall
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Begin CPU_TimeF
//
//
      
	icall = icall+1;
	if(icall == 1)
	{
		cout<<"CPU_TimeF 1.0    9 April     2004' CPU time Uses FORTRAN95 intrinsic timing functions";
		cout<< endl;
	}


//
//
      call SYSTEM_CLOCK(COUNT=itime, COUNT_RATE=irate)
      TIME = dble(itime)/dble(irate)
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  End CPU_TimeF
//
//
//
      RETURN
//
      END SUBROUTINE CPU_TimeF

	  */



void Allocate_Input()
{

	int ii;
	int ierror[5] = {999,999,999,999,999};


	cout<<endl<<"Allocate_Input  1.0   16 May 2014 ,Memory allocation to temporary input arrays"<<endl;

         
	Rk_name = new char[MaxNum_Media]; 
//	MAT=new char[MaxNum_Media];
	MAT=new string[MaxNum_Media];
//	STAT=ierror(40));
	
//cout<< "   1   "<<MAT[0]<<"     2 "<<MAT[1]<<"fin";  We understood that the allocation is done for more than one strings
	
	if (!Rk_name || !MAT)
	{
		ierror[0]=1;
	}
	else
	{
		ierror[0]=0;
	}




	StateI= new int[MaxNum_Media];
	StateIndex = new int[MaxNum_Media];
	
	if (!StateI || !StateIndex)
	{
		ierror[1]=1;
	}
	else
	{
		ierror[1]=0;
	}	





	comm = new char[50];
	if (!comm)
	{
		ierror[2]=1;
	}
	else
	{
		ierror[2]=0;
	}



	XIN = new float [NEQ,MaxNum_Media];
	if (!XIN)
	{
		ierror[3]=1;
	}
	else
	{
		ierror[3]=0;
	}

	YIN= new float [NEQ,MaxNum_Media];
	if (!YIN)
	{
		ierror[4]=1;
	}
	else
	{
		ierror[4]=0;
	}




		 	


	for (ii=0 ; ii<=4 ; ii++)
	{
		if(ierror[ii] == 0)
		{
			cout<<"Memory allocation at point in subroutine Allocate_Input was successful  " << ii+40;
			cout<< endl;
		}
		else if(ierror[ii] == 999)
		{
			continue;
		}
		else
		{
			cout<< "ERROR- Memory allocation at point in subroutine Allocate_Input was unsuccessful";
			cout<<"THE EXECUTION WAS STOPPED    //////////////ERROR"<<ii;
			cout<< endl;

			system("pause");
			exit(0);
		}

	}

	
}



void READ_MainInputF()
{

//***********************************************************************
//***********************************************************************
//*                                                                     *
//*         ROUTINE FOR READING THE DATA IN THE MAIN INPUT FILE         *
//*          DIRECTLY FROM THE "HydrateResSim" INPUT DATA BLOCK         *
//*                                                                     *
//*                  Version 1.0 - September 29, 2004                   *     
//*                                                                     *
//***********************************************************************
//***********************************************************************
//
      //IMPLICIT NONE
// 
// -------
// ... Integer variables
// -------
// 
    //  INTEGER, INTENT(OUT) :: MC


	
	int MC;

// 
    
	int ICALL = 0, icom = 0;
// 
    
	int N_DomInit;
// 
    
	int IE1,IE18,ier1,I_match,ico,i,j,k,n;
// 
// -------
// ... Character variables     
// -------
// 

	string KeyWord;

	char w75[75];
// 
   //   SAVE ICALL,icom
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of READ_MainInputF
//
//

	ICALL = ICALL+1;

	if(ICALL == 1)
	{
		cout<<endl<<"READ_MainInputF V1.0,   21 May 2014, Read all data provided through the main input file"<<endl<<endl;
	}

// 
// 
//***********************************************************************
//*                                                                     *
//*                           INITIALIZATIONS                           *
//*                                                                     *
//***********************************************************************
// 
// 

	DELTMX = 0;  //

	ITI    = 0;     // ITI = # of print-out time data specified
//

	N_DomInit = 0;   // Number of domains for domain-wise initialization
//

	MC = 0;          // Flag for default mesh file (MC /= 0 when the MINC concept is involved)
//

	IRPD = 5  ;      // Default relative permeability (all phases fully mobile)

	ICPD = 1  ;      // Default capillary pressure (PC = 0)
//

	CPD[1] = 0;  // Default capillary pressure parameters (1...3)

	CPD[2] = 0;

	CPD[3] = 1;
//

	iddiag = 0 ;     // Flag for consideration of diffusion (ÒDIFFUÓ-data) 
// 
// ----------
// ... Initialize logical variables        
// ----------



	
	NoFloSimul = false;   // NoFloSimul = .FALSE. : Flow simulation will be ... 
                              // ... conducted - Reset to .TRUE. if flow simulation is bypassed
// 



	NoVersion  = false;  // NoVersion = .FALSE. : Version info is to be ...
                                // ... printed - Reset to .TRUE. if reset by the NOVER keyword
//
      
	SolverBlok = false;   // NoSolvrBlk = .FALSE. : Solver data from MOP(21); 
                                //    Reset to .TRUE. if the SOLVR block is present
//
      
	CoordNeed  = false;     // CoordNeed = .FALSE. : Coordinate arrays are ...  
                                // ... needed - Will be internally reset to .TRUE. if otherwise
//
      
	WellTabData  = false;   // WellTabdata = .FALSE. : There are no tabular well data  
                                //    Will be internally reset to .TRUE. if otherwise 
//
    
	RadiHeat_Flag = false;   // Radiative heat transport is not accounted for 
                                //    Will be internally reset to .TRUE. if otherwise 
//
     
	CondHeat_Flag = true;    // Conductive heat transport is accounted for 
                                //    Will be internally reset to .FALSE. if otherwise 
//
    
	RandmOrdrIC   = false;    // RandmOrdrIC = .FALSE. : Initial conditions are entered ...  
                                 //    in the element order - denoted by the presence of the START keyword 
                                 //    Will be internally reset to .TRUE. (random order) if otherwise 



	// 
// ----------
// ... Initialize character variables        
// ----------
// 
      PermModif = 'NONE';    // PermModif = 'NONE' : No cell-by-cell permeability modification 
// 
// ----------
// ... Initialize primary variable values - CAREFUL// Array operations        
// ----------
// 
	  for (i =1; i<=NK1;i++)
	  {
		  DEP[i] = 0; 
	  }
	  for (i =1; i<=NK1;i++)
	  {
		  XX[i] = 0; 
	  }
//
// ... Initialize fluxes: CAREFUL// Whole array operation
//
      if(NPH > 1)
	  {
		  conV->FluxVD      = 0;


		  for(n=0; n<=MaxNum_Conx-1;n++) 
		  {
			  for (i=0; i<=Num_MobPhs-1;i++)
			  {
				  conV[n].FluxVF[i] =0;
			  }
			  
		  }
	  }
//
      conx->FluxH  = 0;
//
	  for (n=0;n<=MaxNum_Conx-1;n++) 
	  {
		  for (i=0; i<=Num_MobPhs-1;i++)
			  {
				  conx[n].FluxF[i] = 0;
				  conx[n].VelPh[i] = 0;
			  }
	  }

// 
// ----------
// ... Initialize component and equation parameters       
// ----------
// 
      NEQ1  = NEQ+1  ; 
      NK1   = NK+1   ;
      NFLUX = 2*NEQ+1;

// 
//***********************************************************************
//*                                                                     *
//*                 READING AND SEARCHING FOR KEYWORDS                  *
//*                                                                     *
//***********************************************************************
//

GOTOKEYWORD:

	  cin >> KeyWord;

	  if (KeyWord=="ROCK" || KeyWord=="MEDIA")
	  {
		  cout<<endl<<" ROCK KEYWORD WAS ENTERED SUCCESSFULLY! PLEASE IMPORT RELATED VALUES"<<endl;
		  READ_RockData();
		  goto GOTOKEYWORD;
	  }
	  else if (KeyWord=="RPCAP")
	  {
		  cin>> IRPD>> RPD[1]>> RPD[2]>> RPD[3]>> RPD[4]>> RPD[5]>> RPD[6]>> RPD[7]  ;
          cin>> ICPD>> CPD[1]>> CPD[2]>> CPD[3]>> CPD[4]>> CPD[5]>> CPD[6]>> CPD[7]  ;
		  goto GOTOKEYWORD;
	  }
	  else if (KeyWord=="MESHM")
	  {

		  //MESHM(MC);
		  goto GOTOKEYWORD;
	  }
	  else if  (KeyWord=="START")
	  {
         RandmOrdrIC = true;   //Initial conditions to be read in random order
		 cout<<endl<<" START KEYWORD WAS ENTERED SUCCESSFULLY";
		 goto GOTOKEYWORD;
	  } 
	  else if  (KeyWord=="HYDRATE")
	  {
         cout<<endl<<" HYDRATE KEYWORD WAS ENTERED SUCCESSFULLY, PLEASE ENTER RELATED DATA"<<endl<<endl;
		 READ_HydrateData();
		 goto GOTOKEYWORD;
	  }
	  else if  (KeyWord=="PARAM")
	  {
         cout<<endl<<" PARAM KEYWORD WAS ENTERED SUCCESSFULLY, PLEASE ENTER RELATED DATA"<<endl<<endl;
		 cout<< " READ COMPUTATIONAL PARAMETERS"<<endl;
		 READ_ParamData();
		 goto GOTOKEYWORD;
	  }
	  else if  (KeyWord=="ELEMENT")
	  {
         cout<<endl<<" ELEMENT KEYWORD WAS ENTERED SUCCESSFULLY, PLEASE ENTER RELATED DATA"<<endl<<endl;
		 READ_ElemData();
		 goto GOTOKEYWORD;
	  }
	  else if  (KeyWord=="CONNECTION")
	  {
         cout<<endl<<" CONNECTION KEYWORD WAS ENTERED SUCCESSFULLY, PLEASE ENTER RELATED DATA"<<endl<<endl;
		 READ_ConxData();
		 goto GOTOKEYWORD;
	  }
	  else if  (KeyWord=="COFT")
	  {
         
		  cout<<endl<<" COFT KEYWORD WAS ENTERED SUCCESSFULLY, PLEASE ENTER RELATED DATA"<<endl<<endl;
		 
		 //
// ...... Initialization - CAREFUL// Whole array operation
//
         
		  N_obs_conx    = 0;
		  obs_conx_Flag = false;
		  for (i=0;i<=100-1;i++)
		  {
			  obs_conx.name[i]="                ";
			  obs_conx.num[i]  = 0;
		  }
		  
//
// ...... Read in the observation connections
//
         for (i=0;i<=100-1;i++)
		 {
			 cin>>obs_conx.name[i];
             if(obs_conx.name[i] == "ENDOFCOFT")
			 {
				cout<<"COFT FINISHED SUCCESSFULY"<<endl<<endl;
				 break;
			 }
		 }

         N_obs_conx = i-1;
//
         if(N_obs_conx > 0)
		 {
            obs_conx_Flag = true;
            //OPEN(14,file='Conx_Time_Series',status='unknown')
            //REWIND 14
		 }

		 goto GOTOKEYWORD;
	  }
	  else if  (KeyWord=="GENER")
	  {
         cout<<endl<<" GENER KEYWORD WAS ENTERED SUCCESSFULLY!"<<endl<<endl;
		 READ_GnerData();
		 goto GOTOKEYWORD;

	  }
	  else if  (KeyWord=="ENDCI")
	  {
         cout<<endl<<" ENDCI KEYWORD WAS ENTERED SUCCESSFULLY! END OF KEYWORD IMPORT"<<endl<<endl;

	  }
	  else if  (KeyWord=="INIT")
	  {
         cout<<endl<<" INIT KEYWORD WAS ENTERED SUCCESSFULLY! SET THE INITIAL CONDITION!!!"<<endl<<endl;
		 READ_InConData();
		 goto GOTOKEYWORD;
	  }
	  else
	  {

		  // ...... Storing lines with unknown leading keywords 
	  
		  icom = icom+1;


		  if(icom <= 50) 
		  {
			 //comm(icom) = KeyWord //w75
		  }
		  goto GOTOKEYWORD;
	  }
	  
	 
}


void READ_RockData()
{
	/*
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*    ROUTINE FOR READING THE DATA DESCRIBING THE MEDIA PROPERTIES     *
!*         DIRECTLY FROM THE "HydrateResSim" INPUT DATA BLOCK          *
!*                                                                     *
!*                   Version 1.0 - January 6, 2004                     *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
*/

    //  IMPLICIT NONE

// -------
//... Double precision variables
//-------

      float OrgCF,VolGrain;
// 
// -------
// ... Integer variables
// -------
// 
      int ICALL = 0;
      int NAD,n,n_i,i,iiii;
// 
// -------
// ... Character variables
// -------
// 
      string MAT_name;
// 
      //SAVE ICALL
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of READ_RockData
//
//
      ICALL = ICALL+1;
      if (ICALL == 1)
	  {
		   cout<<endl<<"READ_RockData     1.0    22 May   2014, Read the rock type data from the ROCK block of the input data file"<<endl; 
	  }
                                    
//
//
//
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>                                                                     >
//>                         ROCK TYPE LOOP                              >
//>                                                                     >
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
//
//
      //DO_NumRock: DO n_i = 1,MaxNum_Media+1
	  for (n_i = 0;MaxNum_Media;n_i++)
	  {

// 
//***********************************************************************
//*                                                                     *
//*             Checking the adequacy of the dimensioning               *
//*                                                                     *
//***********************************************************************
// 
         if(n_i == MaxNum_Media)
		 {   
			 cin>>MAT_name ;     // Reading the MaxNum_Media+1 record
// 
// .......... If 'MAT_name' is not blank, print an error message & stop 
// 
             if(MAT_name != "EndRK")
			 {
				 cout<<endl<<endl<<" ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR"<<endl;
				 cout<<"----------------------S I M U L A T I O N   A B O R T E D----------------------"<<endl;
				 cout<<"The number of input rock types"<< n_i<<"  exceeds the maximum declared number MaxNum_Media, CORRECT AND TRY AGAIN" ;
				 cout<< endl;
				 system("pause");
				 exit(0);				 

			 }
			 else
			 {
// ......... Otherwise, the rock data read-in is completed 
// 
              
				 N_PFmedia = MaxNum_Media;  // The number of media is set
                 cout<<"ROCK DATA WERE IMPORTED SUCCESSFULLY. END OF ROCK SECTION."<<endl<<endl;
					 break;
			 }
//
		 }

//***********************************************************************
//*                                                                     *
//*              READING THE ROCK TYPE DATA: 1st Record                 *
//*                                                                     *
//***********************************************************************
// 
         PoMed[n_i].RGrain = 0.0;
// 
         cin>> MAT[n_i];
		 cin>> NAD;
		 cin>> media[n_i].DensG;
		 cin>> media[n_i].Poros;
		 cin>> media[n_i].Perm[1];
		 cin>>media[n_i].Perm[2];
		 cin>>media[n_i].Perm[3];
		 cin>> media[n_i].KThrW;
		 cin>> media[n_i].SpcHt;
		 cin>> PoMed[n_i].RGrain;
// 
// --------
// ... End of the rock records        
// --------
// 
		 if(MAT[n_i] == "     ")
		 {
// 
// ......... Determine the number of media in the system 
// 
            N_PFmedia = n_i-1;
// 
// ......... Otherwise, the rock data read-in is completed 
// 
            
			break;
// 
			}


//***********************************************************************
//*                                                                     *
//*             INITIALIZATION OF DATA IN SECOND RECORD                 *
//*                                                                     *
//***********************************************************************
// 
         media[n_i].Compr = 0;  
         media[n_i].Expan = 0;     
         media[n_i].KThrD = 0;    
         media[n_i].Tortu = 0;     
         media[n_i].Klink = 0;     
// 
         media[n_i].NumF_RP = 0;           
         media[n_i].NumF_CP = 0;           
// 
// -----------
// ...... Account for the number of additional records for the cases of  
// ...... (a) reference conditions or (b) seeds for permeability multipliers        
// -----------
//        
         if (MAT[n_i] == "REFCO" || MAT[n_i] == "SEED")
		 {
            if(NAD != 0) 
			{
				NAD = 0;
			}
		 }


// 
//***********************************************************************
//*                                                                     *
//*              READING THE ROCK TYPE DATA: 2nd Record                 *
//*                                                                     *
//***********************************************************************
// 
         if(NAD > 0)
		 {
// 
			 cin>>media[n_i].Compr;
			 cin>> media[n_i].Expan;
			 cin>> media[n_i].KThrD;
			 cin>> media[n_i].Tortu;
			 cin>> media[n_i].Klink;
			 cin>>OrgCF;

// 
		 }
// 
// -----------
// ...... Assign default value to dry thermal conductivity if needed  
// -----------
//        
         if (abs(media[n_i].KThrD) < 1e-6)
		 {
			 media[n_i].KThrD = media[n_i].KThrW;
		 }
// 
//***********************************************************************
//*                                                                     *
//*           READING THE ROCK TYPE DATA: 3rd & 4th Records             *
//*                                                                     *
//***********************************************************************
// 
         if(NAD > 1)
		 {
// 
            cin>>  media[n_i].NumF_RP>> media[n_i].ParRP[1]>> media[n_i].ParRP[2]>> media[n_i].ParRP[3]>> media[n_i].ParRP[4]>> media[n_i].ParRP[5]>> media[n_i].ParRP[6] >> media[n_i].ParRP[7];

            cin>>  media[n_i].NumF_CP>> media[n_i].ParCP[1]>> media[n_i].ParCP[2]>> media[n_i].ParCP[3]>> media[n_i].ParCP[4]>> media[n_i].ParCP[5]>> media[n_i].ParRP[6] >> media[n_i].ParRP[7];
// 
		 }
// 
// ...... Printing the rock type number and name 
// 
         cout<< endl<< "DOMAIN NO."<<n_i<<"     MATERIAL NAME -- "<<MAT[n_i]<<endl;

//
// <<<                      
// <<< End of the ROCK TYPE LOOP         
// <<<
//




	  }
// 



}


void READ_HydrateData()
{


// ...... Modules to be used 
// 
//         USE Basic_Param
// 
//         USE RefRGas_Param
//         USE Hydrate_Param
// 
//         USE RealGas_Properties
//
//***********************************************************************
//***********************************************************************
//*                                                                     *
//*           ROUTINE FOR READING THE HYDRATE PROPERTY DATA             *
//*             DIRECTLY FROM THE TOUGH90 INPUT DATA BLOCK              *
//*                                                                     *
//*                  Version 1.0 - September 30, 2004                   *     
//*                                                                     *
//***********************************************************************
//***********************************************************************
//
//      IMPLICIT NONE
// 
// -------
// ... Double precision arrays
// -------
// 
      float Sum_H;
// 
// -------
// ... Integer variables
// -------
// 
      int ICALL = 0, n,m,ier;
// 
// -------
// ... Character parameters
// -------
// 
      string Reaction_Type;  // Indicates the dissociation type
// 
// -------
// ... Logical parameters
// -------
// 
      bool EX;

	  string EOS_NAME;
// 
      //SAVE ICALL
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of READ_HydrateData
//
//
      ICALL = ICALL+1;
      if(ICALL == 1)
	  {
		  cout<<endl<<"READ_HydrateData  1.0   30 September 2014, Read the hydrate data from the HYDRATE block of the input data file"<<endl;
	  }
// 
//***********************************************************************
//*                                                                     *
//*        OPENING A FILE FOR STORING HYDRATE DISSOCIATION DATA         *
//*                                                                     *
//***********************************************************************
// 

	  ifstream HYD_OUT ("Hydrate_Info");
	  if (HYD_OUT.is_open())
	  {
		   cout<<"FILE 'Hydrate_Info'  EXISTS --- OPEN AS AN OLD FILE"<<endl;
		   HYD_OUT.close();
		   HYD_OUT.open("Hydrate_Info");
	  }
	  else
	  {
		  cout<<"FILE 'Hydrate_Info'  DOES NOT EXIST --- OPEN AS A NEW FILE"<<endl;
		  HYD_OUT.close();
		  ofstream HYD_OUT ("Hydrate_Info");
	  }


// 
//***********************************************************************
//*                                                                     *
//*     READING THE BASIC DATA IDENTIFYING THE HYDRATE COMPOSITION      *
//*                                                                     *
//***********************************************************************
// 
     // cin>> (5,*)HCom.NCom    // The number of hydrate-producing gases
	  cin>> HCom.NCom ;   // The number of hydrate-producing gases
// 
      cout<<endl<<"H Y D R A T E       R E L A T E D       I N F O R M A T I O N"<<endl;         // Write a heading into the output file
// 
//***********************************************************************
//*                                                                     *
//*      INITIALIZE THE MAXIMUM # OF COMPONENTS IN THE GAS PHASE        *
//*                                                                     *
//***********************************************************************
// 
         MaxNum_GasComp = HCom.NCom+1;   // Set equal to the number ...
//                                       // ... of hydrate gases + 1 (H2O) 
//
         GPhase_Com= new string [MaxNum_GasComp];          // Allocate memory to the array 



 

	if (!GPhase_Com)
	{
		ier=1;
	}
	else
	{
		ier=0;
	}


// 
// 
	
	if(ier == 0)
	{
		cout<<"Memory allocation of 'GPhase_Com' in subroutine 'READ_HydrateData' was successful ";
		cout<< endl;
	}
	else
	{
		cout<< "ERROR- Memory allocation of 'GPhase_Com' in subroutine 'READ_HydrateData' was unsuccessful!  THE EXECUTION WAS STOPPED    !!";
		cout<< endl;
		cin.get();
		system("pause");
		exit(0);
	}

	if(HCom.NCom < 1 || HCom.NCom > 2) 
	{
		cout<<endl<<"ERROR -      The number of hydrate-producing gases 'Num_HydrComp' can be either 1'  a simple hydrate or Num_HydrComp = 2, binary hydrates, SIMULATION ABORTED - CORRECT AND TRY AGAIN,ERROR"<< HCom.NCom<<endl;
		system("pause");
		exit(0);
	}

	cout<<endl<<" HYDRATE COMPOSITION, 'Number of Components of the Composite', 'Gas Hydrate = ', 'The components are "<<HCom.NCom<<endl;


	     
	
	  
	
	for (n = 0;n<=HCom.NCom-1;n++)
	{

		cin>> HCom.nameG[n]>> HCom.hydrN[n]>>  HCom.moleF[n] ;  //The mole fraction in the composite hydrate


// ...... Check if available in the data base  
//
		if (HCom.nameG[n] == "CH4" || HCom.nameG[n] == "C2H6" || HCom.nameG[n] == "C3H8" ||  HCom.nameG[n] == "H2S" ||  HCom.nameG[n] == "CO2" ||  HCom.nameG[n] == "N2" )
		{ // Acceptable gas options
			GPhase_Com[n] = HCom.nameG[n] ;           // Make the compound a component in the gas phase
		}
		else
		{
			cout<< endl<< "ERROR-The substance is not among the 'hydrate-producing gases in the ISHR HYDRATE database, SIMULATION ABORTED - CORRECT AND TRY AGAIN',ERROR"<< HCom.nameG[n]<<endl;    //Unavailable gas compound option
			system("pause");
			exit(0);                           //Print a warning and stop
		}

 
	}


//-------- 
// ... Ensure that the mole fractions add up to one  
//-------- 
// 
	
    Sum_H=0;
	for (int iiii=0;iiii<=HCom.NCom-1;iiii++)
	{
	Sum_H += HCom.moleF[iiii] ;            // CAREFUL// Whole array operations
	}

    
	if (abs(1.0e0 - Sum_H) > 1.0e-6)
	{
		cout<< endl<< "ERROR, The sum of the mole fractions of the components of the hydrate  do not add up to one, SIMULATION ABORTED - CORRECT AND TRY AGAIN, ERROR"<< endl;    //Unavailable gas compound option
		system("pause");
		exit(0);                               // Stop the simulation
	}        
//
// ... Make water the last component of the gas phase  
// 
      
	GPhase_Com[HCom.NCom] = "H2O";

// 
//***********************************************************************
//*                                                                     *
//*     Set up the basic parameters for gas behavior description        *
//*                                                                     *
//***********************************************************************
	cout<<endl<<"Equation of State:Peng_Robinson: PR, Soave-Redlick-Knowng:SRK, Redlick-Knowng:RK"<<endl;

    cin>>EOS_NAME; 
	RGasSetup(MaxNum_GasComp, GPhase_Com,EOS_NAME,'G');


//***********************************************************************
//*                                                                     *
//*          Compute the component hydrate molecular weights            *
//*                                                                     *
//***********************************************************************
// 
      for (n=0;n<=HCom.NCom-1;n++)   // The hydrate component loop 
	  {
		  
		  HCom.MolWt[n] =  GasR[ID[n]].MolWt +( GasR[7].MolWt  *HCom.hydrN[n]) ;        // Hydration number !        
//																						// Molecular weight of Water and base gas is needed! 
																						// GasR[7]=Water

		  HCom.GasMF[n] =  GasR[ID[n]].MolWt  / HCom.MolWt[n]  ;          // Molecular weight of component hydrate (kg/mol)

		  HCom.H2OMF[n] =  1.0e0 - HCom.GasMF[n];
// 
// ...... Convert to g/mol   
// 
		  HCom.MolWt[n] = 1.0e3*HCom.MolWt[n];
// 
//----------- 
// ...... Print out the data   
//----------- 
// 
		  cout<<endl<<endl;
		  cout<<"Component #" << n+1 << "-Hydrate                   : "<<   HCom.nameG[n]       <<endl;
		  cout<<"Hydration number                       = "<<HCom.hydrN[n]  <<endl;
		  cout<<"Mole Fraction in the composite hydrate = " << HCom.moleF[n]<<endl;
		  cout<<"Molecular weight of component hydrate  = "<< HCom.MolWt[n] <<" g/mol" <<endl;   // GH-component molecular weight;
// 
// 
	  }

    
	  
//-------- 
// ... Compute the molecular weight of the composite hydrate  
//-------- 
// 


	  MW_Hydrate=0;
	
	  for (int iiii=0;iiii<=HCom.NCom-1;iiii++)
	  {
		  MW_Hydrate += HCom.MolWt[iiii] *HCom.moleF[iiii] ;            // CAREFUL// Whole array operations
	  }	  
	  //MW_Hydrate = SUM( HCom.MolWt(1:HCom%NCom)         *HCom.moleF(1:HCom%NCom)         );

	  cout<<endl<<"PROPERTIES OF THE COMPOSITE HYDRATE"<<endl<<"MOLECULAR WEIGHT                       = "<< MW_Hydrate  <<" g/mol" ;  // Write a heading in the output file  
	  cout<<endl;



//***********************************************************************
//*                                                                     *
//*            READING THE HYDRATE THERMAL CONDUCTIVITY DATA            *
//*                                                                     *
//***********************************************************************
// 
      cin>> GH.N_ThC;  // The number of coefficients in the thermal conductivity ...
                           // ... polynomial F = A0+A1*T+A2*T^2+...+Am*T^m
// 
//-------- 
// ... Allocate memory to the corresponding coefficient array  
//-------- 
// 
      

	  GH.p_ThC = new float[GH.N_ThC];
// 
// ... Print explanatory comments  
// 

	
	  if (!GH.p_ThC)
	  {
		  ier=1;
	  }
	  else
	  {
		  ier=0;
	  }

	  if(ier == 0)
	  {
		  cout<<endl<<"Memory allocation of GH.p_ThC in subroutine READ_HydrateData was successful";
		  cout<< endl;
	  }
	  else
	  {
		  cout<< "ERROR- Memory allocation of GH.p_ThC in subroutine READ_HydrateData was unsuccessful!!!!!    THE EXECUTION WAS STOPPED    !!!!!!! "<<endl;
		  cout<< endl;
		  cin.get();
		  system("pause");
		  exit(0);
	  }


//-------- 
// ... Read the coefficients   
//-------- 
// 
	  
	  for (int iiii=0;iiii<=GH.N_ThC-1;iiii++)
	  {
		  cin>>GH.p_ThC[iiii] ;      // Coefficients
	  }


// 
//-------- 
// ... Print out the data   
//-------- 
// 
      cout<<endl<<endl<<"THERMAL CONDUCTIVITY FUNCTION F_ThC = F(T)"<< endl ;
	  cout<< "Number of coefficients GH.N_ThC in the F(T) = A0+A1*T+A2*T^2+...+Am*T^m polynomial = " << GH.N_ThC <<endl;
	  	  
	  for (int iiii=0;iiii<=GH.N_ThC-1;iiii++)
	  {
		   cout<< "Coefficients A"<<iiii<<" =" << GH.p_ThC[iiii] << "=> Units: W/m/C^"<<iiii+2<<" "<<endl;   // Coefficients       
	  }
	  cout<<endl;
	  
//        
// 


//***********************************************************************
//*                                                                     *
//*              READING THE HYDRATE SPECIFIC HEAT DATA                 *
//*                                                                     *
//***********************************************************************
// 
      cin>> GH.N_SpH  ;    // The number of coefficients in the specific heat ...
                          // ... polynomial F = A0+A1*T+A2*T^2+...+Am*T^m
// 
//-------- 
// ... Allocate memory to the corresponding coefficient array  
//-------- 
// 
	  

	  GH.p_SpH = new float[GH.N_SpH];

	  if (!GH.p_SpH)
	  {
		  ier=1;
	  }
	  else
	  {
		  ier=0;
	  }

// 
// ... Print explanatory comments  
// 

	  if(ier == 0)
	  {
		  cout<<endl<<"Memory allocation of GH.p_SpH in subroutine READ_HydrateData was successful";
		  cout<< endl;
	  }
	  else
	  {
		  cout<< "ERROR- Memory allocation of GH.p_SpH in subroutine READ_HydrateData was unsuccessful!!!!!    THE EXECUTION WAS STOPPED    !!!!!!! "<<endl;
		  cout<< endl;
		  cin.get();
		  system("pause");
		  exit(0);
	  }


// 
//-------- 
// ... Read the coefficients   
//-------- 
//

	  for (int iiii=0;iiii<=GH.N_SpH-1;iiii++)
	  {
		  cin>>GH.p_SpH[iiii] ;      // Coefficients
	  }


// 
//-------- 
// ... Print out the data   
//-------- 
// 
      cout<<endl<<endl<<"SPECIFIC HEAT FUNCTION F_SpH = F(T)"<< endl ;
	  cout<< "Number of coefficients GH.N_SpH in the F(T) = A0+A1*T+A2*T^2+...+Am*T^m polynomial = " << GH.N_SpH <<endl;
	  	  
	  for (int iiii=0;iiii<=GH.N_SpH-1;iiii++)
	  {
		   cout<< "Coefficients A"<<iiii<<" =" << GH.p_SpH[iiii] << "=> Units: J/kg/C^"<<iiii+2<<" "<<endl;   // Coefficients       
	  }
	  cout<<endl;







//***********************************************************************
//*                                                                     *
//*                    READING THE HYDRATE DENSITY                      *
//*                                                                     *
//***********************************************************************
// 
      cin>> GH.N_Rho;  // The number of coefficients in the density ...
                          // ... polynomial F = A0+A1*T+A2*T^2+...+Am*T^m
// 
//-------- 
// ... Allocate memory to the corresponding coefficient array  
//-------- 
//

	  GH.p_Rho = new float[GH.N_Rho];

	  if (!GH.p_Rho)
	  {
		  ier=1;
	  }
	  else
	  {
		  ier=0;
	  }


// 
// ... Print explanatory comments  
// 
 	  if(ier == 0)
	  {
		  cout<<endl<<"Memory allocation of GH.p_Rho in subroutine READ_HydrateData was successful";
		  cout<< endl;
	  }
	  else
	  {
		  cout<< "ERROR- Memory allocation of GH.p_Rho in subroutine READ_HydrateData was unsuccessful!!!!!    THE EXECUTION WAS STOPPED    !!!!!!! "<<endl;
		  cout<< endl;
		  cin.get();
		  system("pause");
		  exit(0);
	  }
// 
//-------- 
// ... Read the coefficients   
//-------- 
// 

	  for (int iiii=0;iiii<=GH.N_Rho-1;iiii++)
	  {
		  cin>>GH.p_Rho[iiii] ;      // Coefficients
	  }


// 
//-------- 
// ... Print out the data   
//-------- 
//
	        
	  cout<<endl<<endl<<"DENSITY FUNCTION F_Rho  = F(T)"<< endl ;
	  cout<< "Number of coefficients GH.N_Rho in the F(T) = A0+A1*T+A2*T^2+...+Am*T^m polynomial = " << GH.N_Rho <<endl;
	  	  
	  for (int iiii=0;iiii<=GH.N_Rho-1;iiii++)
	  {
		   cout<< "Coefficients A"<<iiii<<" =" << GH.p_Rho[iiii] << "=> Units: kg/m^3/C^"<<iiii+1<<" "<<endl;   // Coefficients       
	  }
	  cout<<endl;
 // 
//***********************************************************************
//*                                                                     *
//*     READING THE TEMPERATURE SHIFTS DUE TO SALT AND INHIBITORS       *
//*                                                                     *
//***********************************************************************
// 
// ... Inhibitor-related temperature shift in the hydrate equilibrium curve
// 
      cin>> T_MaxOff >>  C_MaxOff >>   MW_Inhib >>  D_Inhib >>  H_InhSol >>   DifCo_Inh ; // Diffusion coefficient of inhibitor (m^2/s)
// 
//-------- 
// ... Print out the data   
//-------- 
// 
	  cout<< endl;
      cout<<"INHIBITOR-INDUCED TEMPERATURE SHIFT IN THE EQUILIBRIUM HYDRATION CURVE (Salt or Alcohol)"<<endl;
	  cout<<"Maximum temperature shift T_MaxOff                         = " << T_MaxOff <<" (degrees C)"<<endl;
	  cout<<"Mole fraction C_MaxOff at which the maximum T-shift occurs = " << C_MaxOff << endl;
	  cout<<"Molecular weight of the inhibitor MW_Inhib                 = " << MW_Inhib<< " (g/mol)" << endl;
	  cout<<"Density of the inhibitor D_Inhib                           = " << D_Inhib << " (kg/m^3)" << endl;
	  cout<<"Enthalpy of inhibitor dissolution D_InhSol                 = " << H_InhSol<< " (J/kg)" << endl;
	  cout<<"Diffusion coefficient DifCo_Inh of inhibitor in H2O        = " << DifCo_Inh << "(m^2/s)" << endl;


//***********************************************************************
//*                                                                     *
//*            READING THE HYDRATE EQUATION OPTION NUMBER               *
//*                                                                     *
//***********************************************************************
// 
      cin>> F_EqOption ;  // Flag for selecting the hydrate P-T and H-T functions
//                            // = 0: Moridis [2003]; = 1: Kamath (1984)
//-------- 
// ... Print out the data   
//-------- 
// 
      cout<< endl << "HYDRATION EQUATIONS OPTION NUMBER, F_EqOption = " << F_EqOption << " Flag for selecting the equations used for (a) the hydration P-T and (b) H-T relationhips     = 1: The Kamath [1984] equations are used =1: Extended forms of the Moridis [2003] equations are used (DEFAULT)" ;
// 
// 
//***********************************************************************
//*                                                                     *
//*              READING THE TYPE OF HYDRATE DISSOCIATION               *
//*                                                                     *
//***********************************************************************
// 
      
	  cout<<endl;
	  cin>> Reaction_Type;
// 
// ... For kinetic dissociation 
// 
      if(Reaction_Type == "KINETIC")
	  {
		  F_Kinetic = true;
		  cout<<endl<< "TYPE OF HYDRATE FORMATION OR DISSOCIATION : K I N E T I C , PARAMETERS OF KINETIC DISSOCIATION"<<endl;
	  }
	  else if (Reaction_Type == "EQUILIBRIUM")
	  {
// 
// ... For equilibrium dissociation 
// 
		  F_Kinetic = false;
          cout<<endl<< "TYPE OF HYDRATE FORMATION OR DISSOCIATION : E Q U I L I B R I U M " << endl;
	  }
	  else
	  {
// ... On error, print a message and stop 
// 
         
		  cout << "ERROR - The variable Reaction_Type is neither EQUILIBRIUM or KINETIC "<< endl;
		  cout << "Unable to determine type of hydrate reaction, SIMULATION ABORTED - CORRECT AND TRY AGAIN -ERROR" ;
		  cout<< endl;
		  cin.get();
		  system("pause");
		  exit(0);
	  }


//----------- 
// ...... Ensure correct dimensioning   
//----------- 
// 
      if(F_Kinetic == true)
	  {       
         if(NK < 3) 
		 {
            cout<<endl;
			cout<<"ERROR- The number of mass components NK = " << NK << " is in conflict with the kinetic type of hydrate reaction. For a kinetic reaction, NK must be either 3 (without salt) or 4 (with salt)"<<endl;
			cout<<"SIMULATION ABORTED - CORRECT AND TRY AGAIN, ERROR"<<endl;
			cin.get();
			system("pause");
			exit(0);
		 }

//         
         if(NEQ < 4)
		 {
            cout<<endl;
			cout<<"ERROR, The number of equations NEQ = " << NEQ << " is in conflict with the kinetic type of hydrate reaction. For a kinetic reaction, NEQ must be either 4 (without salt) or 5 (with salt) " << endl;
			cout<< "SIMULATION ABORTED - CORRECT AND TRY AGAIN 'ERROR" <<endl; 
			cin.get();
			system("pause");
			exit(0);
		 }

	  }




//***********************************************************************
//*                                                                     *
//*                     FOR KINETIC HYDRATION ONLY                      *
//*                                                                     *
//***********************************************************************
// 
      if(F_Kinetic == true)
	  {
// 
//----------- 
// ...... Read and write the kinetic parameters   
//----------- 
// 
        
		  cin>> Activ_En  >>  InRate_k0 >>  Area_Factor  ;     // Adjustment factor for reaction area = ratio of ...
//                                  // ... actual grain surface area/spherical surface area
         
		  cout<< "Activation energy       = " << Activ_En<< " J/mol" <<endl;
		  cout<< "Intrinsic rate constant = " <<  InRate_k0<< " mol/(m^2.Pa.s)" << endl;
		  cout<< "Area adjustment factor  = " <<  Area_Factor << endl;


// 
// ...... Convert the intrinsic dissociation energy to [kg/(m^2.Pa.s)] 
// 
//         InRate_k0 = 1.0d3*InRate_k0/MW_Hydrate
         
		  InRate_k0 = InRate_k0*MW_Hydrate/1.0e3;
// 
//----------- 
// ...... Allocate memory to arrays related to kinetic dissociation  
//----------- 
// 
         SArea_Hyd = new float [MaxNum_Elem];
		 
// 
// ...... Print explanatory comments  
// 

	  
		 if (!SArea_Hyd)
		 {
			 ier=1;
		 }
		 else
		 {
			 ier=0;
		 }

		 if(ier == 0)
		 {
			 cout<<endl<<"Memory allocation of SArea_Hyd in subroutine READ_HydrateData was successful";
			 cout<< endl;
		 }
		 else
		 {
			 cout<< "Memory allocation of SArea_Hyd in subroutine READ_HydrateData was unsuccessful!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!! ERROR"<<endl;
			 cout<< endl;
			 cin.get();
			 system("pause");
			 exit(0);
		 }
// 
// 
// 
	  }
// 
//***********************************************************************
//*                                                                     *
//*            Allocate memory to double precision array of             *
//*         inhibitor mass fraction at source (injection well)          *
//*                                                                     *
//***********************************************************************
// 
// 
      Inh_WMasF = new float [MaxNum_SS];

// 
// ... Print explanatory comments  
// 
	  if (!Inh_WMasF)
	  {
		  ier=1;
	  }
	  else
	  {
		  ier=0;
	  }

	  if(ier == 0)
	  {
		  cout<<endl<<"Memory allocation of Inh_WMasF in subroutine READ_HydrateData was successful";
		  cout<< endl;
	  }
	  else
	  {
		  cout<< "ERROR, Memory allocation of Inh_WMasF in subroutine READ_HydrateData was unsuccessful!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!ERROR"<<endl;
		  cout<< endl;
		  cin.get();
		  system("pause");
		  exit(0);
	  }


// 
// ... Print dividers for improved legibility     
// 

	  cout<< "*******************************************************************************"<<endl;
	  cout<< "*******************************************************************************"<<endl;

}

















void RGasSetup(int N_Com,string Cname[],string EquTyp,char PhaseF)
{
 
//         USE RefRGas_Param
//
//***********************************************************************
//***********************************************************************
//*                                                                     *
//*         Providing reference to the built-in gas parameters          *     
//*                                                                     *
//*                    Version 1.00 - 8 April, 2014                     *     
//*                                                                     *
//***********************************************************************
//***********************************************************************
//
//      IMPLICIT NONE
// 
// -------
// ... Integer variables
// -------
// 
//      INTEGER, INTENT(IN) :: N_Com
// 
      int ICALL = 0;
// 
      int i,j;
// 
// -------
// ... Character arrays
// -------
// 
//      char   Cname[6][N_Com];
// 
// -------
// ... Character variables
// -------
// 
//      char EquTyp[3];
//      char PhaseF;
// 
      //SAVE ICALL
// 
// 
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of RGasSetup
// 
// 
      ICALL = ICALL+1;
	  if(ICALL == 1)
	  {
		  cout<< endl<< "RGasSetup         1.0    8 April     2014', 'Mapping to the built-in reference ', 'parameter data set for real gas behavior"<<endl;
	  }

// 
// -------
// ... Define the module-specific parameters
// -------
// 
      NumCom   = N_Com;      
      EOS_Type = EquTyp;
      Phase    = PhaseF;
// 
// -------
// ... Check if requested phase is available
// -------
// 
	  if (Phase != 'L' && Phase != 'G')   // If unavailable phase, ...
	  {
		  cout<< " R U N   A B O R T E D, The phase indicator in parameter PHASE  "<< Phase<<"   must be either L (Liquid), or G (Gas), 'CORRECT AND TRY AGAIN', ERROR"<< endl;   // ... write a message, and ... 
		  system("pause");
		  exit(0);                           //Print a warning and stop                                           // ... abort the simulation
	  }
// 
//***********************************************************************      
//*                                                                     *
//*                       Select the EOS to use                         *
//*                                                                     *
//***********************************************************************      
//
      if(EOS_Type=="SRK")
	  {   // Soave-Redlich-Kwong
       //  CALL K_SRK
	  }
	  else if (EOS_Type=="PR")
	  {        // Peng-Robinson
       //  CALL K_PR
	  }
	  else if(EOS_Type=="RK")
	  {    // Redlich-Kwong
       //  CALL K_RK
	  }
	  else
	  {// If unavailable EOS, ...
         cout<<" The type of EOS equation EOS_Type = "<< EOS_Type <<" is not an available option  !!!! ====>   THE RUN WAS ABORTED - CORRECT AND TRY AGAIN')" <<endl;
         // ... write a message, and ...
		 system("pause");
		 exit(0);                           //Print a warning and stop 
	  }
// 
//***********************************************************************      
//*                                                                     *
//*          Determine the indices for communication with the           *
//*                    built-in parameter data base                     *
//*                                                                     *
//***********************************************************************      
//
      iH2O = 0;  // Initializations
      iH2S = 0;
// 
// 
// 
      
	  
	  for (j=0;j<=NumCom-1;j++)
	  {
// 
// ----------
// ...... Determining the ID indices       
// ----------
// 
       
		  
		  if (Cname[j] == "CH4")
		  {
			  ID[j] = 0;
		  }
		  else if (Cname[j] == "C2H6")
		  {
            ID[j] = 1;
		  }
		  else if (Cname[j] == "C3H8")
		  {
			  ID[j] = 2;
		  }
		  else if (Cname[j] == "H2S")
		  {
			  ID[j] = 3;
			  iH2S  = 1;
		  }
		  else if (Cname[j] == "CO2")
		  {
			  ID[j] = 4;
		  }
		  else if (Cname[j] == "N2")
		  {
			  ID[j] = 5;
		  }
		  else if (Cname[j] == "O2")
		  {
			  ID[j] = 6;
		  }
		  else if (Cname[j] == "H2O")
		  {
			  ID[j] = 7;
			  iH2O  = 1;
		  }
		  else if(Cname[j] == "C2H3OH")
		  {
            ID[j] = 8;
		  }
		  else
		  {
            cout<<" The component "<< Cname[j]<<"   is not an available option !!!!! ==>   THE RUN WAS ABORTED - CORRECT AND TRY AGAIN'"<<endl; 
         // ... write a message, and ...
		    system("pause");
		 	exit(0);                           //Print a warning and stop 
		  }
// 
// ----------
// ...... Mapping the mixture parameters to the reference ones      
// ----------
// 
         Gas[j].TCrit = GasR[ID[j]].TCrit;        
         Gas[j].PCrit = GasR[ID[j]].PCrit;        
         Gas[j].VCrit = GasR[ID[j]].VCrit;       
         Gas[j].ZCrit = GasR[ID[j]].ZCrit; 
//      
         Gas[j].Omega = GasR[ID[j]].Omega;       
         Gas[j].MolWt = GasR[ID[j]].MolWt;       
         Gas[j].DMom  = GasR[ID[j]].DMom ;      
//      
	  }
// 
//***********************************************************************      
//*                                                                     *
//*         Mapping the EOS coefficients to the reference ones          *
//*                                                                     *
//***********************************************************************      
//
      for (i=0;i<=NumCom-1; i++)
	  {
		  for (j=0;j<=NumCom-1; j++)
		  {
			  if (EOS_Type=="SRK")
			  {
				  Kij_SC[i][j] = Kij_S[ID[i]][ID[j]];
			  }
			  else if (EOS_Type == "PR")
			  {
				  Kij_PC[i][j] = Kij_P[ID[i]][ID[j]];
			  }
			  else
			  {
				  Kij_RC[i][j] = Kij_R[ID[i]][ID[j]];
			  }
		  }
	  }

			   


//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of RGasSetup
// 
//
      //     
	  
} //END of RGasSetup
//


void READ_ParamData()
{

// 
// ...... Modules to be used 
// 
//         USE Basic_Param
//         USE GenControl_Param
//         USE Diffusion_Param
// 
//         USE Variable_Arrays
//
//***********************************************************************
//***********************************************************************
//*                                                                     *
//*       ROUTINE FOR READING THE MAIN COMPUTATIONAL PARAMETERS         *
//*                                                                     *
//*                 Version 1.0 - September 21, 2004                    *     
//*                                                                     *
//***********************************************************************
//***********************************************************************
//
//      IMPLICIT NONE
// 
// -------
// ... Integer variables
// -------
// 
      int ICALL = 0;
      int I,J,N;
// 
// -------
// ... Character variables
// -------
// 
      //char State_id[3];
	  string State_id;
// 
//    SAVE ICALL
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of READ_ParamData
//
//
      ICALL = ICALL+1 ;
      if(ICALL == 1)
	  {
//		  WRITE(11,6000);
		  cout<< "READ_ParamData    1.0   3 June 2014', Read the computational parameter data from the PARAM block of the input data file"<<endl;
	  }
	  	  
	  
	  //WRITE(6,6050)
	  //cout<<"READ Parameters Data!!! PARAM block of the input data file In 3 Steps"<<endl;


// 
//***********************************************************************
//*                                                                     *
//*       READING THE COMPUTATIONAL PARAMETER DATA: 1st Record          *
//*                                                                     *
//***********************************************************************
// 
      
	  cin>> NOITE  >>   KDATA >>  MCYC   >>   MSEC  >>  MCYPR;
	  
	  for (int iiii=0;iiii<=23;iiii++)
	  {
		  cin>>    MOP[iiii];
	  }
	  cin >> DIFF0  >>  TEXP >>  BE;
//      
// ... Reset MOP(14) to 1 (T-dependent function) for hydrate simulation  
// 
      if(EOS_Name == "HYDRATE" || MOP[13]!= 1)
	  {
		  MOP[13] = 1 ;
	  }

// 
// ... Allow more than 10000 time steps 
// 
      
	  if(MCYC < 0)
	  {
		  MCYC = 100*abs(MCYC);
	  }

// 
//***********************************************************************
//*                                                                     *
//*       READING THE COMPUTATIONAL PARAMETER DATA: 2nd Record          *
//*                                                                     *
//***********************************************************************
// 
      cin>> TSTART >> TIMAX >> DELTEN >> DELTMX >> ELST >> GF >> REDLT >> SCALE;
// --------
// ... For DELTEN < 0, read list of specific DELTs to be used         
// --------
// 
      NDLT = 0;
// 
      if(DELTEN < 0.0e0)
	  {
		  NDLT = -(DELTEN);
		  for (int N=0; N<=NDLT-1;N++)
		  {
			  for (int J=0;J<=7;J++)
			  {
				  cin >> DLT[J+8*N];
			  }
		  }

// 
          

		  DELTEN = DLT[1];
// 
	  }
//

//***********************************************************************
//*                                                                     *
//*       READING THE COMPUTATIONAL PARAMETER DATA: 3rd Record          *
//*                                                                     *
//***********************************************************************
// 
     cin>> RE1 >>  RE2 >>  U_p >>  WUP >>   WNR >>  DFAC >>  FOR >> State_id ;
// ----------
// ... Determine the state index
// ----------
// 
      if (State_id=="Gas")
	  {
// 
// ... Single phase: Gas           
//    
      
		  GenStIndx = 1;
	  }
	  else if  (State_id=="Aqu")
	  {
// ... Single phase: Aqueous           
		  GenStIndx = 2;
	  }
	  else if  (State_id=="AqG")
	  {
//...Two phases: Aqueous+Gas   
//    
		  GenStIndx = 3;
	  }
	  else if (State_id=="IcG")
	  {
// 
// ... Two phases: Ice+Gas   

         GenStIndx = 4 ;
	  }
	  else if (State_id=="AqH")
	  {
// 
// ... Two phases: Aqueous+Hydrate    
		  GenStIndx = 5;
	  }
	  else if (State_id=="IcH")
	  {
// 
// ... Two phases: Ice+Hydrate    
         GenStIndx = 6;
	  }
	  else if (State_id=="AGH")
	  {
// 
// ... Three phases: Aqueous+Hydrate+Gas    
//    
         GenStIndx = 7;
	  }
	  else if (State_id=="AIG")
	  {

// ... Three phases: Aqueous+Ice+Gas    
         GenStIndx = 8;
	  }
	  else if (State_id=="AIH")
	  {
// ... Three phases: Aqueous+Ice+Hydrate    
         GenStIndx = 9;
	  }
	  else if (State_id=="IGH")
	  {
// ... Three phases: Gas+Ice+Hydrate    
         GenStIndx = 10;
	  }
	  else if (State_id=="QuP")
	  {

// 
// ... Quadruple point: Gas+Ice+Hydrate+Aqueous    
// 
         GenStIndx = 11;
	  }
	  else
	  {

// 
// ... Otherwise: Error ...    
         cout<< "ERROR-                      S I M U L A T I O N   A B O R T E D         The state identifier of uniform initial conditions State_id_U = " <<  State_id << "This is not an available option for this Equation of State. CORRECT AND TRY AGAIN!!!! ERROR" << endl;
		 system("pause");
		 exit(0); 
	  }


	  // --------
// ... Assign default values        
// --------
// 
      if(NOITE == 0) NOITE = 8;
      if(KDATA == 0) KDATA = 1;
      if(MCYPR == 0) MCYPR = 1;
// 
//      IF(TEXP  == 0.0d0) TEXP  = 1.8d0
      if(REDLT == 0.0e0) REDLT = 4.0e0;
      if(SCALE == 0.0e0) SCALE = 1.0e0;
// 
      if(RE1 == 0.0e0) RE1 = 1.0e-5;
      if(RE2 == 0.0e0) RE2 = 1.0e0;
      if(FOR == 0.0e0) FOR = 1.0e0;
// 
      if(WUP <= 0.0e0) WUP = 1.0e0;
      if(WNR == 0.0e0) WNR = 1.0e0;
// 
      if(U_p  == 0.0e0) U_p  = 1.0e-1;
// 
// --------
// ... Initialize some parameters        
// --------
// 
      KCYC  = 0;        // Number of elapsed time steps
      ITER  = 0;        // Number of Newtonian iterations in current Dt
      ITERC = 0;        // Cumulative number of Newtonian iterations
      TIMIN = TSTART;   
// 
//***********************************************************************
//*                                                                     *
//*        READING THE INITIAL VALUES OF THE PRIMARY VARIABLES          *
//*                                                                     *
//***********************************************************************
// 
      for (int i=0;i<=NK1-1;i++)
	  {
		  cin>>DEP[i];
	  }


}

void READ_ElemData()
{
//	     SUBROUTINE READ_ElemData
// 
// ...... Modules to be used 
// 
//         USE Basic_Param
// 
//         USE PFMedProp_Arrays, ONLY: MAT
//
//***********************************************************************
//***********************************************************************
//*                                                                     *
//*     ROUTINE FOR READING THE DATA DESCRIBING THE DOMAIN ELEMENTS     *
//*         DIRECTLY FROM THE "HydrateResSim" INPUT DATA BLOCK          *
//*                                                                     *
//*                   Version 1.0 - January 6, 2004                     *     
//*                                                                     *
//***********************************************************************
//***********************************************************************
//
//      IMPLICIT NONE
// 
// -------
// ... Double precision variables
// -------
// 
      float VOLX,AHTX,zref,X,Y,Z;
// 
// -------
// ... Integer variables
// -------
// 
      int ICALL = 0;
// 
      int NE,NSEQ,NSEQ1,NADD,MA2;
      int n,i,m,i_ASCII;
// 
// -------
// ... Character variables
// -------
// 
//      char EL[3],MA1[3];
	  string EL,MA1;
      //char MA12[5];
      //char EL6_0[6];
	   
	  string MA12;
      string EL6_0;
// 
// -------
// ... Logical variables
// -------
// 
      bool Medium_ByName;
// 
      //SAVE ICALL
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of READ_ElemData
//
//
      ICALL = ICALL+1     ;
      if(ICALL == 1)
	  {
		  cout<<"READ_ElemData     1.0    4 June   2014', Read the element data from the ELEME block of the input data file";
		  cout<<endl;

	  }

	  //WRITE(6,6000)
// 
// --------
// ... Printing headings        
// --------
// 
      //PRINT 5002
      //REWIND 4
// 
// --------
// ... Printing according to the number of characters in the element names         
// --------
// 
      if(No_CEN != 8)
	  {
//
// ................ 5-Character elements
//
                   
		  cout << "WRITE FILE *MESH* FROM INPUT DATA"<<endl;
                   //WRITE(4,6002)
	  }
	  else
	  {
       
// ................ 8-Character elements
		  cout<< "WRITE FILE *MESH* FROM INPUT DATA -',        ' Mesh data follow the ',a4,' option" << " ext8"<<endl;
       //            WRITE(4,6006) 'ext8'
	  }
// 
//***********************************************************************
//*                                                                     *
//*                       READING THE ELEMENT DATA                      *
//*                                                                     *
//***********************************************************************
// 
GOTOELEM:
	 
	  if(No_CEN != 8)
	  {
//
// ............... Read data for 5-Character elements
//
		  cin>> EL >>  NE >>  NSEQ  >>   NADD >>  MA12 >>  VOLX >>   AHTX >>  zref >>  X >>  Y >>  Z ;
	  }
	  else
	  {
//
// ............... Read data for 8-Character elements
//
		  cin>> EL6_0   >>   NE >>  MA12 >>    VOLX >>  AHTX >>   zref >>  X >>  Y >>  Z ;
//         
		  //EL   = EL6_0(1:3);

		  EL[0]   = EL6_0[0];
		  EL[1]   = EL6_0[1];
		  EL[2]   = EL6_0[2];

		  NSEQ = 0;
		  NADD = 0;
	  }
// 
// --------
// ... End of the element records         
// --------
// 
	  if(EL=="ELD")
	  {
		  cout<<"END OF ELEMENT INPUT";
		  return;

	  }
//
      NSEQ1 = NSEQ+1;
// 
//***********************************************************************
//*                                                                     *
//*                  ASSIGNING THE MATERIAL/ROCK TYPE                   *
//*           WRITING THE ELEMENT DATA INTO THE ÒMESHÓ FILE             *
//*                                                                     *
//***********************************************************************
// 
      Medium_ByName = false;                     // Default: Medium by number
// 
      for (int i=0;i<=4;i++)
	  {
		  //i_ASCII = int(MA12[i:1]);               // Determine the ASCII value of ...
//                                                 // ... each character of the medium name
		  i_ASCII = int(MA12[i]);               
		  if(i_ASCII == 32)
		  {
			  continue;       // Continue iteration for blanks
//                                                 
		  }

		  if(i_ASCII > 57 || i_ASCII <48)
		  {   // If any character is a letter, ...
			  Medium_ByName = true;                // ... the flag is reset
			  break;
		  }
 
	  }

// 
// -----------
// ...... Elements specified by material number         
// -----------
// 
      if(Medium_ByName = false)
	  {
// 
       
		  cin>> MA12>> MA2;        // Read MA2(target) from MA12(string) 
//
// ...... 5-Character elements
//
         
		  if(No_CEN != 8)
		  {
            for (i=0;i<=NSEQ1-1;i++)
			{
				n = NE+(i-1)*NADD;

               	cout<< EL<<n<<MA2<<VOLX<<AHTX<<zref<<X<<Y<<Z;
			}
		  }
		  else
		  {
//
// ...... 8-Character elements
//
             for (i=0;i<=NSEQ1-1;i++)
			 {
				 n = NE+(i-1)*NADD;
                 cout<<EL6_0 << n << MA2 << VOLX << AHTX << zref << X << Y <<Z;
			 }
		  }
//
// ...... Continue reading data
//
	  }

      goto GOTOELEM;
}

void READ_ConxData()
{

 
// ...... Modules to be used 
// 
//         USE Basic_Param
//   
//         USE Connection_Arrays, ONLY: conx
//
//***********************************************************************
//***********************************************************************
//*                                                                     *
//*   ROUTINE FOR READING THE DATA DESCRIBING THE DOMAIN CONNECTIONS    *
//*         DIRECTLY FROM THE "HydrateResSim" INPUT DATA BLOCK          *
//*                                                                     *
//*                   Version 1.0 - January 6, 2004                     *     
//*                                                                     *
//***********************************************************************
//***********************************************************************
//
//      IMPLICIT NONE
// 
// -------
// ... Double precision variables
// -------
// 
      float D1,D2,AREAX,BETAX,sigx;
// 
// -------
// ... Integer variables
// -------
// 
      int ICALL = 0;
// 
      int NE1,NE2,NSEQ,NSEQ1,NAD1,NAD2,ISOT;
      int N,N1,N2,I;
// 
// -------
// ... Character variables
// -------
// 
//      char EL1[3],EL2[3];

	  string EL1,EL2;

//      char EL6_1[6],EL6_2[6];
	  
	  string EL6_1,EL6_2;
// 
      //SAVE ICALL
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of READ_ConxData
//
//
      ICALL = ICALL+1;
      if(ICALL == 1)
	  {
		  cout<< "READ_ConxData     1.0    18 June   2014, Read the connection data from the CONNE block of the input data file";
	  }

// 
// --------
// ... Printing headings and initializing        
// --------
// 

	  
	  //      WRITE(4,6002)    //%%%%%%%%%%%%%%%%%%%%%%%  WRITE TO THE FILE

	  
	  //
// ... Initializing the connection counter to zero
//
      NCON  = 0;
// 
//***********************************************************************
//*                                                                     *
//*                     READING THE CONNECTION DATA                     *
//*                                                                     *
//***********************************************************************
//

//***********************************************************************
//*                                                                     *
//*                     READING THE CONNECTION DATA                     *
//*                                                                     *
//***********************************************************************
// 
GOTOCONNECTION:


	  if(No_CEN != 8)
	  {
//
// ............... Read connection data for 5-Character elements
//
       
		  cin>> EL1 >> NE1 >> EL2 >> NE2 >> NSEQ >> NAD1 >> NAD2 >> ISOT >> D1 >> D2 >> AREAX >> BETAX >> sigx ;
	  }
	  else
	  {
//
// ............... Read connection data for 8-Character elements
//
       
		  
		  cin>> EL6_1 >>  NE1 >> EL6_2 >>  NE2 >> ISOT >> D1 >> D2 >>  AREAX >> BETAX >> sigx;

		  EL1[0]  = EL6_1[0];
		  EL1[1]  = EL6_1[1];
		  EL1[2]  = EL6_1[2];

               
		  NSEQ = 0;
		  NAD1 = 0;
		  NAD2 = 0;

	  }
// 
// --------
// ... End of the connection records - No info on connection elements        
// --------
// 
      
	  
	  
	  if(EL1 == "ECL" && NE1 == 0)
	  {

// 
		  // ................. The connection data read-in is completed 
// 
           return;
// 
	  }

	  
// --------
// ... End of the connection records - Info on connection elements        
// --------
// 
      if(EL1 == "+++")
	  {
// 
       //             WRITE(4,5002)
// 
// ................. Read the connection-associated element numbers 
// 
                    
		  if(No_CEN != 8)
		  {
//cin>>  (conx(N)%n1,conx(N)%n2,N=1,NCON);
			  
			  for (N=0;N<=NCON-1;N++)
			  {
				  cin>> conx[N].n1 >> conx[N].n2;
//WRITE(4,5010) (conx(N)%n1,conx(N)%n2,N=1,NCON)

			  }
		  }
		  else
		  {
//        READ 5015,    (conx(N)%n1,conx(N)%n2,N=1,NCON)
			  
			  for (N=0;N<=NCON-1;N++)
			  {
				  cin>> conx[N].n1 >> conx[N].n2;
//WRITE(4,5010) (conx(N)%n1,conx(N)%n2,N=1,NCON)

			  }

//                      WRITE(4,5015) (conx(N)%n1,conx(N)%n2,N=1,NCON)
		  }
// 
// ................. The connection data read-in is completed 
// 
       //             WRITE(4,5001)
                    
		  return;
// 
	  }
//
      NSEQ1 = NSEQ+1;
// 
//***********************************************************************
//*                                                                     *
//*          WRITING THE CONNECTION DATA INTO THE ÒMESHÓ FILE           *
//*                                                                     *
//***********************************************************************
// 
      for (int i=0;i<=NSEQ1-1;i++)
	  {
//
// ...... Count connections
//
         
		  NCON = NCON+1;
//
// ...... Determine element numbers in connection NCON
//
         
		  N1 = NE1+(i-1)*NAD1;
         
		  N2 = NE2+(i-1)*NAD2;
//
// ...... Writing the element info into MESH
//
         
		  if(No_CEN != 8)
		  {

//			  WRITE(4,6004) EL1,N1,EL2,N2,ISOT,D1,D2,                    AREAX,BETAX,sigx
		  }
		  else
		  {
          
//			  WRITE(4,6006) EL6_1,NE1,EL6_2,NE2,ISOT,D1,D2,                    AREAX,BETAX,sigx
		  }
		  // 
	  
	  
	  }
// 
// --------
// ... Continue reading connection data
// --------
// 
      goto GOTOCONNECTION;








}

void READ_GnerData()
{


// ...... Modules to be used 
// 
         //USE Basic_Param
// 
         //USE Q_Arrays
//
//***********************************************************************
//***********************************************************************
//*                                                                     *
//*     ROUTINE FOR READING THE GENERATION (SOURCE AND SINK) DATA       *
//*         DIRECTLY FROM THE "HydrateResSim" INPUT DATA BLOCK          *
//*                                                                     *
//*                  Version 1.0 - September 11, 2004                   *     
//*                                                                     *
//***********************************************************************
//***********************************************************************
//
      //IMPLICIT NONE
// 
// -------
// ... Double precision variables
// -------
// 
      float GX,EX,HX,Inh_WInjMF;
// 
// -------
// ... Integer variables
// -------
// 
      int ICALL = 0;
// 
      int NE,NS,N1,N2,NSEQ,NSEQ1,NADD,NADS,LTAB;
      int n,i,k,LTABA,ier;
// 
// -------
// ... Character variables
// -------
// 
      char ITAB;
      string EL,SL;
      string SS_TypX;
      string EL6_0;
// 
      //SAVE ICALL
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of READ_GnerData
//
//
      ICALL = ICALL+1   ;
      if(ICALL == 1) 
	  {
		  cout<<"READ_GnerData     1.0   19 June 2014, Read the source/sink data from the GENER block of the input data file";
		  cout<<endl<<endl;
	  }

// 
// --------
// ... Printing headings and attaching GENER file        
// --------
// 
//      REWIND 3
//
      cout<<"WRITE FILE *GENER* FROM INPUT DATA"<<endl;
 
	 
	  ofstream GENER_OUT;
	  GENER_OUT.open("GENER", ios::app);
	  GENER_OUT << "GENER\n";
	  GENER_OUT << "START TO WRITE IN THE 'GENER' FILE\n";
	  GENER_OUT.close();
//
// ... Seeting the source/sink counter to zero
//
      NOGN  = 0;

// 
//***********************************************************************
//*                                                                     *
//*                    READING THE SOURCE/SINK DATA                     *
//*                                                                     *
//***********************************************************************
// 
GOTOGENER:

	  if(No_CEN != 8)
	  {
//
// ............... Read data for 5-Character elements
//
           
		  cin>> EL >> NE >> SL >> NS>> NSEQ >> NADD >> NADS >> LTAB >> SS_TypX >> ITAB >> GX >> EX >> HX >> Inh_WInjMF; 
	  }
	  else
	  {
//
// ............... Read data for 8-Character elements
//
          
		  cin>> EL6_0 >> NE >> SL >>NS >> NSEQ >> NADD >> NADS >> LTAB >> SS_TypX >> ITAB >> GX >> EX >> HX>>Inh_WInjMF;

		  EL[0]   = EL6_0[0];
		  EL[1]   = EL6_0[1];
		  EL[2]   = EL6_0[2];
          NSEQ = 0;
          NADD = 0;
          NADS = 0;
	  
	  }
// 
// --------
// ... End of the source/sink records - No info on correlation to elements        
// --------
// 
      if(EL == "EOG" && NE == 0)
	  {
		  GENER_OUT << "      ";
		  
// 
// ................. The source/sink data read-in is completed 
// 
          
		  return;
	  }

	  
// 
// --------
// ... End of the source/sink records - Info on correlation to elements        
// --------
// 
      if(EL == "+++")
	  {
		  
		  GENER_OUT << "+++   ";
	
// 
// ................. Read the source/sink-associated element numbers 
// 
    
		  if(No_CEN != 8)
		  {

			  for (i=0;i<=NOGN-1;i++)
			  {
				  cin>> SS[n].el_num;
			  }
			  for (i=0;i<=NOGN-1;i++)
			  {
				  GENER_OUT<< SS[n].el_num;
			  }
		  }
		  else
		  {
			  
			  for (i=0;i<=NOGN-1;i++)
			  {
				  cin>> SS[n].el_num;
			  }

			  for (i=0;i<=NOGN-1;i++)
			  {
				  GENER_OUT<< SS[n].el_num;
			  }
		  }
// 
// ................. The source/sink data read-in is completed 
// 
                    
		  GENER_OUT<<"                 ";
          return;
//
					
	  }


     
	  NSEQ1 = NSEQ+1;
      
	  if(LTAB == 0)
	  {
		  LTAB=1;
	  }
// 
//***********************************************************************
//*                                                                     *
//*         WRITING THE CONNECTION DATA INTO THE 'GENER' FILE           *
//*                                                                     *
//***********************************************************************
// 
      for(i=0;i<=NSEQ1-1;i++)
	  {

//
// ...... Count sources/sinks
//
         
		  NOGN = NOGN+1;
//
// ...... Determine number component of element names 
//
         
		  N1 = NE+(i-1)*NADD;
          N2 = NS+(i-1)*NADS;
//
// ...... Writing the element info into file GENER
//
         if(No_CEN != 8)
		 {
			 GENER_OUT<< EL<< N1<< SL<< N2<< LTAB<< SS_TypX<< ITAB<<  GX<< EX<< HX<< Inh_WInjMF;
		 }
		 else
		 {
            GENER_OUT<<EL6_0<< N1<< SL<< N2<< LTAB<< SS_TypX<< ITAB<< GX<< EX<< HX<< Inh_WInjMF;
		 }
// 
        
		 LTABA = abs(LTAB);
// 
//***********************************************************************
//*                                                                     *
//* READING TABULAR DATA FROM THE 'GENER' DATA BLOCK IN THE INPUT FILE  *
//*                                                                     *
//***********************************************************************
// 
         if(LTABA > 1 && SS_TypX != "DELV")
		 {
// 
// -------------
// ......... Allocate memory to the temporary arrays F1, F2, F3    
// -------------
// 
          
			 F1=new float [LTABA];
			 F2=new float [LTABA];
			 F3=new float [LTABA];

// ......... Print explanatory comments  
// 
	
			 if (!F1 || !F2 || !F3)
			 {
				 ier=1;
			 }
			 else
			 {
				 ier=0;
			 }

			 

			 if(ier == 0)
	  		 {
				 cout<<endl<<"Memory allocation to arrays F1,F2,F3 in subroutine READ_GnerData was successful";
				 cout<< endl;
			 }
			 else
			 {
				 cout<< "ERROR, Memory deallocation to arrays F1,F2,F3 in subroutine READ_GnerData was unsuccessful!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!ERROR"<<endl;
				 cout<< endl;
				 cin.get();
				 system("pause");
				 exit(0);
			 }

// 
// -------------
// ......... Reading the data    
// -------------
// 
			 if(i< 2)
			 {
				 for (k=0;k<=LTABA-1;k++)
				 {
					 cin>> F1[k] ;  // Tabular data of times
				 }
				 for (k=0;k<=LTABA-1;k++)
				 {
					 cin>> F2[k] ;  // Tabular data of injection/production rates
				 }
              
              
				 if(ITAB != 42) // asci(*)=42
				 {
					 	
					 for (k=0;k<=LTABA-1;k++)
					 {
						 cin>> F3[k] ;  // Tabular data of corresponding specific enthalpy
					 }
					
				 }
			 }

// 
//***********************************************************************
//*                                                                     *
//*             WRITING TABULAR DATA INTO THE 'GENER' FILE              *
//*                                                                     *
//***********************************************************************
// 
			 		
			 for (k=0;k<=LTABA-1;k++)
			 {			 
				 GENER_OUT << F1[k] ; 
			 }

			 for (k=0;k<=LTABA-1;k++)
			 {			 
				 GENER_OUT << F2[k] ;  
			 }

          
 				 
			 if(ITAB != 42)  // asci(*)=42
			 {
				 	
				 for (k=0;k<=LTABA-1;k++)
				 {
					 GENER_OUT << F3[k] ;
				 }
					
			 }

			 
// 
// -------------
// ......... Deallocate memory from the temporary arrays F1, F2, F3    
// -------------
// 
            delete F1;
			delete F2;
			delete F3;
			

// 
// ......... Print explanatory comments  
//
			
			if (!F1 && !F2 && !F3)
			{
				ier=0;
			}
			else
			{
				ier=1;
			}

			
			if(ier == 0)
			{
				cout<<endl<<"Memory allocation to arrays F1,F2,F3 in subroutine READ_GnerData was successful";
				cout<< endl;
			}
			else
			{
				cout<< "ERROR, Memory deallocation to arrays F1,F2,F3 in subroutine READ_GnerData was unsuccessful!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!ERROR"<<endl;
				cout<< endl;
				cin.get();
				system("pause");
				exit(0);
			}
// 
		 }
//
// <<<                      
// <<<           
// <<<
//
	  }
// 
// -----------
// ...... Continue reading source/sink data
// -----------
// 
      goto GOTOGENER;
			  


}

void READ_InConData()
{

//         USE Basic_Param
// 
//         USE Variable_Arrays,   ONLY: XX
//
//***********************************************************************
//***********************************************************************
//*                                                                     *
//*          ROUTINE FOR READING THE INITIAL CONDITION DATA             *
//*        DIRECTLY FROM THE "HydrateResSim" INPUT DATA BLOCK           *
//*                                                                     *
//*                 Version 1.0 - September 26, 2004                    *     
//*                                                                     *
//***********************************************************************
//***********************************************************************
//
//      IMPLICIT NONE
// 
// -------
// ... Double precision variables
// -------
// 
      float PORX,TSTX,TIMINX;
// 
// -------
// ... Integer variables
// -------
// 
      
	  int ICALL = 0;
// 
      int NE,N1,NSEQ,NSEQ1,NADD,KCYCX,ITERCX,NMX;
      int I,J;
// 
// -------
// ... Character variables
// -------
// 
      string EL,State_id;
      string EL6_0;
// 
      //SAVE ICALL
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of READ_InConData
//
//
      ICALL = ICALL+1;
      if(ICALL == 1)
	  {
	
		  cout<<"READ_InConData    1.0   20 June 2014, Read the initial condition data from the INIT block of the input data file";

	  }
// 
// --------
// ... Printing headings and attaching the INCON file        
// --------
// 
//      PRINT 5001
//      PRINT 6002
//
      cout<<"WRITE FILE *INCON* FROM INPUT DATA"<<endl;
 
	 
	  ofstream INCON_OUT;
	  INCON_OUT.open("INCON", ios::app);
	  INCON_OUT << "INITIAL CONDITION\n";
	  INCON_OUT << "START TO WRITE IN THE 'INCON' FILE\n";
	  INCON_OUT.close();
// 
//***********************************************************************
//*                                                                     *
//*            READING THE INITIAL CONDITION DATA: 1st record           *
//*                                                                     *
//***********************************************************************
// 

GOTOINIT:

	  if(No_CEN != 8)
	  {
//
// ............... Read initial condition data for 5-Character elements
//
       cin>> EL >> NE >> NSEQ >> NADD >> PORX >> State_id;
	  }
	  else
	  {
//
// ............... Read initial condition data for 8-Character elements
//
       
		  cin>> EL6_0 >> NE >> NSEQ >> NADD >> PORX >> State_id;
//
       
		  EL[0]   = EL6_0[0];
		  EL[1]   = EL6_0[1];
		  EL[2]   = EL6_0[2];
		  NSEQ = 0;
		  NADD = 0;
	  }

// 
// --------
// ... End of the initial condition records - No continuation info         
// --------
// 
      
	  if(EL == "EIC" && NE == 0)
	  {
		  INCON_OUT<< "               ";
	  
// 
// ................. The initial condition data read-in is completed 
// 
      return;
	  }
// 
       
// --------
// ... End of the initial condition records - Info on continuation        
// --------
// 
      if(EL == "+++")
	  {
		   INCON_OUT<< "++++";
// 
// ................. Read the continuation data 
// 
                   

		  cin>> KCYCX >> ITERCX >> NMX >> TSTX >> TIMINX;

// 
// ................. Write the continuation data into file INCON
// 
          
		  INCON_OUT<< KCYCX << ITERCX << NMX << TSTX << TIMINX;
// 
// ................. The initial condition data read-in is completed 
// 
                    
		  
		  INCON_OUT<< "               ";
          
			  return;
		  
// 
	  }
// 
//***********************************************************************
//*                                                                     *
//* READING THE INITIAL CONDITION DATA: 2nd record (Primary Variables)  *
//*                                                                     *
//***********************************************************************
// 
      for (I=0; I<= NK1-1; I++)
	  {
		  cin>>  XX[I];
	  }

// 
      NSEQ1 = NSEQ+1;
// 
//***********************************************************************
//*                                                                     *
//*      WRITING THE INITIAL CONDITION DATA INTO THE ÒINCONÓ FILE       *
//*                                                                     *
//***********************************************************************
// 
	  for (I=0;I<=NSEQ1-1;I++)
	  {
//
// ...... Determine number component of element name 
//
         N1 = NE+(I-1)*NADD;
//
// ...... Writing the element initial condition info into file INCON
//
         if(No_CEN != 8)
		 {
			 INCON_OUT<< EL<< N1<< PORX<< State_id;
		 }
		 else
		 {

			 INCON_OUT<< EL6_0<< N1<< PORX<< State_id;
		 }
// 
		 for (J=0;J<=NK1-1;J++)
		 {
			 INCON_OUT<< XX[J];
		 }
// 
	  }
// 
// -----------
// ...... Continue reading initial condition data
// -----------
// 
      goto GOTOINIT;
				 
//




}

void NumFloPo_Digits()
{

// 
// ...... Modules to be used 
// 
//         USE GenControl_Param
//
//***********************************************************************
////**********************************************************************
//*                                                                     *
//*    ROUTINE FOR CALCULATING THE NUMBER OF SIGNIFICANT DIGITS FOR     *
//*    FLOATING POINT PROCESSING; ASSIGN DEFAULT FOR DFAC, AND PRINT    *
//*      APPROPRIATE WARNING WHEN MACHINE ACCURACY IS INSUFFICIENT      *
//*                                                                     *
//*                   Version 1.0 - March 11, 2003                      *     
//*                                                                     *
//***********************************************************************
//***********************************************************************
//
//      IMPLICIT NONE
// 
// -------
// ... Double precision variables
// -------
// 
      float DF;
// 
// -------
// ... Integer variables
// -------
// 
      int N10;
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> Main body of NumFloPo_Digits
//
//
      cout<<"NumFloPo_Digits   1.0   11 March     2003, Calculate number of significant digits for floating point arithmetic"<<endl<<endl;
// 
// -------
// ... Machine accuracy from FORTRAN95 intrinsic functions          
// -------
// 
      //N10 = streamsize precision() const;//(1.0e-12)
	  N10=8;
      DF  = sqrt(numeric_limits<double>::epsilon( ));
//
      cout<<"EVALUATE FLOATING POINT ARITHMETIC. FLOATING POINT PROCESSOR HAS APPROXIMATELY " <<N10 << " SIGNIFICANT DIGITS. DEFAULT VALUE OF INCREMENT FACTOR FOR NUMERICAL DERIVATIVES IS DFAC = " << DF << endl;
// 
// -------
// ... Determine the derivative increment         
// -------
// 
      if(DFAC == 0.0e0) 
	  {
		  DFAC = DF;
          cout<<" DEFAULT VALUE FOR DFAC WILL BE USED"<< endl;
	  }
	  else
	  {
		  cout << "USER-SPECIFIED VALUE DFAC = " << DFAC << " WILL BE USED"<< endl;
	  }


// 
// -------
// ... Additional print-outs         
// -------
// 
      if(N10 <= 12 && N10 > 8) 
	  {
		  cout<<"WWWWWWWWWW WARNING WWWWWWWWWW: NUMBER OF SIGNIFICANT DIGITS IS MARGINAL; EXPECT DETERIORATED CONVERGENCE BEHAVIOR"<<endl;
	  }
	  
	  if(N10 <= 8)
	  {
		  cout<< "WWWWWWWWWW WARNING WWWWWWWWWW: NUMBER OF SIGNIFICANT DIGITS IS INSUFFICIENT; CONVERGENCE WILL BE POOR OR FAIL; WWWWWWWWWWWWWWWWWWWWWWWWWWWWW: CODE SHOULD BE RUN IN DOUBLE PRECISION!"<<endl;
	  }


}


void READ_Files()
{










}
