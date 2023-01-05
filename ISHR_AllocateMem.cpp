//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>                                                                     >
//>                                                                     >
//>           ISHR_AllocateMem.cpp: Code unit including the             >
//>         routines allocating memory to all ISHR simulations          >
//>                                                                     >
//>                                                                     >
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//

//
//***********************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***********************************************************************
//


#include <iostream>
#include <fstream>
#include <string>
using namespace std;


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
	int ierror = 99999;

// 
// ----------
// ...... Integer variables     
// ----------
//
	int NREDM_t,MPRIM_t,ii;
	int i,j;
// 
// ----------
// ...... Character variables     
// ----------
// 
	//char B_name = "               ";
	string B_name = "               ";
	char HEADR[26];
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
	cout<<"ASHR_Allocate_MemGen   1.0   1 APR 2014 Allocate memory to most arrays based on size provided by input parameters";
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
    if (HEADR == "ISHR MEMORY ALLOCATION")
	{
	}
	else
	{
		cout<<"ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-  '       S I M U L A T I O N   A B O R T E D',/,     T20, The header 'ISHR MEMORY ALLOCATION' is missing at the top of the input file"<<endl<<"CORRECT AND TRY AGAIN',";
		exit(0);
	}


//
  
	cin >> B_name;
    //EOS_Name = TRIM(ADJUSTL(B_name));
	string EOS_Name = B_name;

//
//
//***********************************************************************
//*                                                                     *
//*  READ THE NUMBER OF COMPONENTS, EQUATIONS, PHASES, SEC. VARIABLES   *
//*                                                                     *
//***********************************************************************
		 

 // 
//

	cin << NK,NEQ,NPH,M_BinDif,M_add;
//
// ...... Number of supplementary parameters 
//

	if(M_add <= 0)
	{
		M_Add = 1             
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
	
	EOS_DefaultNum(ADJUSTR(EOS_Name),  MaxNum_MComp,MaxNum_Equat,   MaxNum_Phase,MaxNum_MobPh,M_Add);

// ...... Determine the number of mobile phases 
//
         
	Num_MobPhs = MaxNum_MobPh;
	NPH        = MaxNum_Phase;

		 // ...... For hydrates, ensure that NEQ = NK+1 
//

	if(EOS_Name(1:7) == "HYDRATE") && (NEQ /= (NK+1))
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

	if (NK  > MaxNum_MComp) ||  (NEQ > MaxNum_Equat) ||   (NPH > MaxNum_Phase)
	{
		cout<<"//,20('ERROR-'),//,T33,  '       S I M U L A T I O N   A B O R T E D',//,     T10,'One or more of the parameters NK,NEQ,NPH = (',  i2,',',i2,',',i2,')',/,     T10,'exceed the maximum values ',      '"MaxNum_MComp", "MaxNum_Equat", ',       '"MaxNum_Phase" = (',i2,',',i2,',',i2,')',/,  /,T32,'               CORRECT AND TRY AGAIN',        //,20('ERROR-'))";
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

	cin>> MaxNum_Elem,MaxNum_Conx,No_CEN,FLG_con;
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
		ierror[1]=1;
	}
	else
	{
		ierror[1]=0;
	}


// 
// ----------
// ...... Allocate memory to double precision allocatable arrays     
// ----------
// 

	AI= new float[MaxNum_Elem];    // ALLOCATE(AI(MaxNum_Elem),      STAT=ierror(2))
	if (!AI)
	{
		ierror[2]=1;
	}
	else
	{
		ierror[2]=0;
	}


	AHT= new float[MaxNum_Elem];    // ALLOCATE(AHT(MaxNum_Elem),     STAT=ierror(3))
	if (!AHT)
	{
		ierror[3]=1;
	}
	else
	{
		ierror[3]=0;
	}
// 

	Pm= new float[MaxNum_Elem];    // ALLOCATE(Pm(MaxNum_Elem),      STAT=ierror(4))
	if (!Pm)
	{
		ierror[4]=1;
	}
	else
	{
		ierror[4]=0;
	}


	X_Coord= new float[MaxNum_Elem];    // ALLOCATE(X_Coord(MaxNum_Elem),    Y_Coord(MaxNum_Elem),    Z_Coord(MaxNum_Elem), STAT=ierror(5))
	Y_Coord= new float[MaxNum_Elem];
	Z_Coord= new float[MaxNum_Elem];
		 
	if (!X_Coord) || (!Y_Coord) || (!Z_Coord)
	{
		feerror[5]=1;
	}
	else
	{
		ierror[5]=0;
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
	if (!conx)
	{
		ierror[6]=1;
	}
	else
	{
		ierror[6]=0;
	}

// 

	if(MaxNum_Phase > 1) 
	{
		conV= new VaporFlux[MaxNum_Conx]; //ALLOCATE(conV(MaxNum_Conx), STAT=ierror(7))
	}
	if (!conV)
	{
		ierror[7]=1;
	}
	else
	{
		ierror[7]=0;
	}
		 // 
// ----------
// ...... Allocate memory to arrays within the derived type: Conx  
// ----------
// 
	for (i = 1; i<=MaxNum_Conx;1)
	{
			 //ALLOCATE(Conx(i).FluxF(Num_MobPhs),    Conx(i).VelPh(Num_MobPhs),     STAT = ierror(8))
		Conx(i).FluxF= new float[Num_MobPhs];
		Conx(i).VelPh= new float[Num_MobPhs];
	}


	if (!Conx(i).FluxF) || (!Conx(i).VelPh)
	{
		ierror[8]=1;
	}
	else
	{
		ierror[8]=0;
	}
// 
// ----------
// ...... Allocate memory to double precision allocatable arrays     
// ----------
// 

	sig= new float[MaxNum_Conx];   //ALLOCATE(sig(MaxNum_Conx), STAT=ierror(9))
	if (!sig) 
	{
		ierror[9]=1;
	}
	else
	{
		ierror[9]=0;
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
		ierror[10]=1;
	}
	else
	{
		ierror[10]=0;
	}


	DX= new float [MPRIM_t];           //ALLOCATE(DX(MPRIM_t),          STAT=ierror(11))
	if (!DX)
	{
		ierror[11]=1;
		
	}
	else
	{
		ierror[11]=0;
	}


	DELX= new float [MPRIM_t];         //ALLOCATE(DELX(MPRIM_t),        STAT=ierror(12))
	if (!DELX)
	{
		ierror[12]=1;
	}
	else
	{
		ierror[12]=0;
	}
//

	R= new float [NEQ*MaxNum_Elem+1];  //ALLOCATE(R(NEQ*MaxNum_Elem+1), STAT=ierror(13))
	if (!R)
	{
		ierror[13]=1;
	}
	else
	{
		ierror[13]=0;
	}
// 

	DOLD= new float [NREDM_t];         //ALLOCATE(DOLD(NREDM_t),        STAT=ierror(14))
	if (!DOLD)
	{
		ierror[14]=1;
	}
	else
	{
		ierror[14]=0;
	}
// 
// ----------

// 
// ----------
// ...... Allocate memory to derived-type arrays - Secondary variables  
// ----------
// 

	Cell_V=new SecParam[MaxNum_Elem,0:NEQ]//ALLOCATE(Cell_V(MaxNum_Elem,0:NEQ), STAT=ierror(15))
	if (!Cell_V)
	{
		ierror[15]=1;
	}
	else
	{
		ierror[15]=0;
	}
// 
// ----------
// ...... Allocate memory for arrays within the derived type  
// ----------
// 

	for (j = 0; j<=NEQ;1)
	{
		for (i=1; i<= MaxNum_Elem;1)
		{
			Cell_V(i,j).p_Satr= new float[MaxNum_Phase];
			Cell_V(i,j).p_KRel = new float [MaxNum_Phase];
			Cell_V(i,j).p_Visc = new float [MaxNum_Phase];
			Cell_V(i,j).p_Dens = new float [MaxNum_Phase];
			Cell_V(i,j).p_Enth = new float [MaxNum_Phase];
			Cell_V(i,j).p_PCap = new float [MaxNum_Phase]
			Cell_V(i,j).p_MasF = new float [MaxNum_MComp,MaxNum_Phase];
				 
			if (!Cell_V(i,j).p_Satr) || (!Cell_V(i,j).p_KRel) || (! Cell_V(i,j).p_Visc) || (!Cell_V(i,j).p_Dens)|| (!Cell_V(i,j).p_Enth) || (! Cell_V(i,j).p_PCap) || (!Cell_V(i,j).p_MasF)
			{
				ierror[16]=1;
			}
			else
			{
				ierror[16]=0;
			}
		}
	}



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

	if(M_BinDif != 0) && (MaxNum_MComp >= 2)
	{
		for (j = 0; j<=NEQ;1)
		{
			for (i = 0; i<=MaxNum_Elem;1)
			{
				Cell_V(i,j).p_DifA = new float[MaxNum_Phase];
				Cell_V(i,j).p_DifB = new float[MaxNum_Phase];  //ALLOCATE(Cell_V(i,j)%p_DifA(MaxNum_Phase),     Cell_V(i,j)%p_DifB(MaxNum_Phase),    STAT = ierror(17));

				if (!Cell_V(i,j).p_DifA) || (!Cell_V(i,j).p_DifB)
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

	if(M_add > 0)
	{
		for (j = 0; j<=NEQ;1)
		{
			for (i=1; i<=MaxNum_Elem;1)
			{
				Cell_V(i,j).p_AddD = new float [M_add];

				if (!Cell_V(i,j).p_AddD)
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
		ierror[19]=1;
	}
	else
	{
		ierror[19]=0;
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
		ierror[20]=1;
	}
	else
	{
		ierror[20]=0;
	}
// 

	if(MaxNum_Phase > 1) 
	{
		for (1 = 1; i<=MaxNum_SS;1)
		{
				 //ALLOCATE(SS(i)%frac_flo(Num_MobPhs),    SS(i)%rate_phs(Num_MobPhs),  STAT = ierror(21));

			SS(i).frac_flo = new float [Num_MobPhs];
			SS(i).rate_phs = new float [Num_MobPhs];

			if (!SS(i).frac_flo) || (!SS(i).rate_phs)
			{
				ierror[21]=1;
			}
			else
			{
				ierror[21]=0;
			}
		}
		 }
// 

	Well_TData= new Tabular[MaxNum_SS];//ALLOCATE(Well_TData(MaxNum_SS), STAT = ierror(21))
	if (!Well_TData)
	{
		ierror[21]=1;
	}
	else
	{
		ierror[21]=0;
	}
// 

	WDel = new Tab_BHP [MaxNum_SS];//ALLOCATE(WDel(MaxNum_SS),       STAT = ierror(22))
	if (!WDel)
	{
		ierror[22]=1;
	}
	else
	{
		ierror[22]=0;
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
		ierror[23]=1;
	}
	else
	{
		ierror[23]=0;
	}
// 

	media= new PFMedium [MaxNum_Media];//ALLOCATE(media(MaxNum_Media), STAT=ierror(24))

	if (!media)
	{
		ierror[24]=1;
	}
	else
	{
		ierror[24]=0;
	}


// 
// 
	for (ii=0;ii<=24;i++)
	{
		if(ierror(ii) == 0)
		{
			cout<<" Memory allocation at point ',i3,  ' in subroutine "Allocate_MemGen" was successful')" << ii
		}
		else if(ierror(ii) == 99999_
		{
			continue
		}
		else
		{
			cout<< "('ERROR-'),//,   T2,' Memory allocation at point ',i3,  ' in subroutine "Allocate_MemGen" was unsuccessful',//,    T2,'!!!!!!!    THE EXECUTION WAS STOPPED    !!!!!!!',  //,20('ERROR-'))"<<ii;
			exit(0);
			
		}
		
	}

	 

}

