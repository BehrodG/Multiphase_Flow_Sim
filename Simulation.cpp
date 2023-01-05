#include <iostream>
#include <fstream>
#include <string>
using namespace std;

extern int KC;

extern int IGOOD;
// 

extern float SUMTIM;
extern float TIMIN;

extern float TIMOUT;

extern float DELT;
extern float DELTEN;
// 

extern int NST;
extern int NK1;
extern int ITER;
extern int IHLV;

extern bool Converge_Flag;
extern bool ConvFail_Flag;

void EOS();

void SIMULATION_Cycle()
{
	// ... Modules to be used 
// 
//      USE Basic_Param
//      USE GenControl_Param
// 
//      USE Element_Arrays
//      USE Variable_Arrays
//
//***********************************************************************
//***********************************************************************
//*                                                                     *
//*          EXECUTIVE ROUTINE THAT PERFORMS THE SIMULATION             * 
//*                     WHILE ADVANCING IN TIME                         *
//*                                                                     *
//*                  Version 1.0 - January 9, 2004                      *     
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

float PNIT,TCUR,Dt_solv,T_dd;

// 
// -------
// ... Integer variables
// -------
// 
int ICALL = 0;

// 
int n,NSTL,KIT;
int n_t,n_t0,i_n,NIT,NIT1;
// 
// -------
// ... Logical variables
// -------
// 
bool ConvFail_Flag,DtCut_Flag,PrintNow_Flag;
	  
	  
	  //SAVE ICALL
//
//
//  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>  Main body of SIMULATION_Cycle
//
//
ICALL = ICALL+1;

if(ICALL == 1) 
{
	//WRITE(11,6000);
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
  
KC    = 0;
IGOOD = 0;
// 

SUMTIM = TIMIN;

TIMOUT = SUMTIM;

DELT   = DELTEN;
// 

NSTL = (NST-1)*NK1;
// 

ITER = 0;

IHLV = 0;

KIT  = 0;
// 




Converge_Flag = false ;   // Flag denoting convergence of the Newton-Raphson ...
                                 // ... iteration - Simulation continuing 
// 

ConvFail_Flag = false ;   // Flag denoting convergence failure ...
                                 // ... leading to abrupt simulation end
// 
// ----------
// ... Initialize thermophysical properties        
// ----------
// 
EOS();          // Calling the Equation-of-State routine



}


void EOS()
{

	system("pause");

}