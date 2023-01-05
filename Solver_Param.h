//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
//
#ifndef Solver_Param 
#define Solver_Param


// 
         //SAVE
// 
// ----------
// ...... Double precision parameters
// ----------
// 
         
float   seed = 1.0e-25;
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
// 
//
#endif

