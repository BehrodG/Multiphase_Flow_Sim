//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//
//
//
#ifndef Q_Arrays 
#define Q_Arrays

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
   
	float  frac_flo;
    
	float  rate_phs;
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

float    F1,F2,F3;
// 
// ----------
// ...... Derived-type arrays     
// ----------
// 

SourceSink  SS;
// 
 
Tabular     Well_TData;
// 
 
Tab_BHP     WDel;
// 
// 
#endif

