//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
//
#ifndef GenControl_Param 
#define GenControl_Param
// 
         //SAVE
// 
// ----------
// ...... Double precision arrays     
// ----------
// 

float DLT, TIS;
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

int   INUM,IPRintMCYC,MCYPR,MSEC;

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

char  PermModif[4];

char  ELST[8];

// 
// ----------
// ...... Common logical variables     
// ----------
// 

bool   CoordNeed,WellTabData,NoFloSimul,NoVersion;

bool   PrintOut_Flag,Converge_Flag,RandmOrdrIC;

bool   RadiHeat_Flag,CondHeat_Flag;
// 
#endif


