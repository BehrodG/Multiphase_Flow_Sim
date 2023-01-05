//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
//
#ifndef TimeSeriesP_Param 
#define TimeSeriesP_Param

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
 
	char name[16][100] ; // Names of observation connections

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
 
Observation_SS     obs_SS ;    // The observation SS list
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


#endif
