//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
///
#ifndef Solver_Param 
#define Solver_Param


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

float    *AI,*AHT,*pm;

float    *X_coord;

float    *Y_coord;

float    *Z_coord;
// 
// ----------
// ...... Derived-type arrays     
// ----------
// 
 
Element  *elem;
// 
// 
#endif

