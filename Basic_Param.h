#ifndef Basic_Param 
#define Basic_Param

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

char GPhase_Com[6]; // Names of the gas components  
// 
// ----------
// ...... Common char variables
// ----------
// 

char FLG_con[6];

char EOS_Name[15];

char TITLE[80];
// 
// ----------
// ...... Common bool variables
// ----------
// 

bool   EX_FileVERS;
// 
// 
#endif

