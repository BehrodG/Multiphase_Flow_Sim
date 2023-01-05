#ifndef EOS_DefaultParam 
#define EOS_DefaultParam
//#include <string>
//#include <iostream>

//char EOS_state[3][11] = ( 'Gas', 'Aqu', 'AqG',  'IcG', 'AqH', 'IcH',  'AGH', 'AIG', 'AIH', 'IGH', 'QuP' );
//char EOS_state[3][11] = { {"Gas"}, {"Aqu"}, {"AqG"},  {"IcG"}, {"AqH"}, {"IcH"},  {"AGH"}, {"AIG"}, {"AIH"}, {"IGH"}, {"QuP"} );
	
string EOS_state[11]={"Gas", "Aqu", "AqG",  "IcG", "AqH", "IcH",  "AGH", "AIG", "AIH", "IGH", "QuP" };


//          CONTAINS


void EOS_DefaultNum (EOS_name,MaxNum_MComp,MaxNum_Equat, MaxNum_Phase,MaxNum_MobPh,M_Add)
{

int MaxNum_MComp,MaxNum_Equat,  MaxNum_Phase,MaxNum_MobPh,M_Add;

// -------------
// ...... Character variables
// -------------

char EOS_name[15];
//
//
// =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of EqTemp_HydrCH4 
//
//

cout<< "EOS_DefaultNum    1.0   29 September 2004',6X, 'Default parameters of the hydrate equation of state";
//
//
//***********************************************************************
//*                                                                     *
//*         DEFAULT NUMBERS FOR THE HYDRATE EQUATION OF STATE           *
//*                                                                     *
//***********************************************************************
// 
// 
       if(EOS_Name(1:11) == "HYDRATE-EQU")
	   {
		   MaxNum_MComp = 3;   // Maximum number of mass components
           MaxNum_Equat = 4;   // Maximum number of equations
           MaxNum_Phase = 4;   // Maximum number of phases
           MaxNum_MobPh = 2;   // Maximum number of mobile phases
	   }
	   else if(EOS_Name(1:11) == "HYDRATE-KIN")
	   {
		   MaxNum_MComp = 4;  // Maximum number of mass components
           MaxNum_Equat = 5;  // Maximum number of equations
           MaxNum_Phase = 4;  // Maximum number of phases
           MaxNum_MobPh = 2;  // Maximum number of mobile phases
	   }
	   else
	   {
		   cout<<"       S I M U L A T I O N   A B O R T E D',/, T20,'The equation of state specified by EOS_Name = ',  A15,' is unavailable'  CORRECT AND TRY AGAIN )";
		   exit(0);
	   }

	   //
// ----------
// ...... Adjust the parameter controlling the size of the auxilliary parameter array
// ----------
//
         if((EOS_Name(1:7) == "HYDRATE") && (M_Add < 3)) 
		 {
			 M_Add = 3;
		 }

}





#endif