#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "libPP_5.0.h"
#include "EntSalArb_MP.h"
#include "GNA.h"

main(){	
///////////////////////////Inicializa parametros de la simulacion
int NDX=100;
float fracON = 1.0;
int NDY=NDX;
int T_max = 400;
int NoEnsambles=110*4;

float Birth1= 2.0;
float Coagulation1= 0.00002; //Brown usa: 0.00002; 
float CoaIntra= 0.0008; //Modelo J-C 0.0008
float Dead1= 1.0;
int RadioBirth1= 10;
int RadioCoa1= 5;
int RadioCoaIntra1= 5;  //Modelo Heteromyopia 20

int CantidadEspecies=50;

float Birth2= 2.0;
float Coagulation2= 0; //Brown usa: 0.00002; 
float Dead2= 1.0;
int RadioBirth2= 10;
int RadioCoa2= 5;



SetBirth(3.0,0);
SetCoagulation(1.0,0);
SetDead(1.0,0);
SetCoagulationIntra(0.0,0);
//SetRadioBirth(RadioBirth1,NoEspecie);
//SetRadioCoa(RadioCoa1,NoEspecie);
EscalaTiempoMetabolico(0);

int CantEspecies;
	for(CantEspecies=1;CantEspecies<=CantidadEspecies;CantEspecies++)
	{
			SetSpecie2(CantEspecies, Birth1, Coagulation1, CoaIntra, Dead1, RadioBirth1, RadioCoa1, RadioCoaIntra1);
	}

//	SetSpecie(1, Birth1, Coagulation1, Dead1, RadioBirth1, RadioCoa1);
//	SetSpecie(2, Birth2, Coagulation2, Dead2, RadioBirth2, RadioCoa2);

//omp_set_num_threads(4);
////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:

		
//Dist_MP MP_RhoDist_1;
//float TamParticion=0.0001;
//float Xini=0.0;
//float Xfin=0.12;
	//InicializaDist_MP(&MP_RhoDist_1, TamParticion,Xini,Xfin);
	
Float1D_MP MP_Correlacion_1;
	InicializaFloat1D_MP(&MP_Correlacion_1, NDX);
			
Float2D_MP MP_Corr2D_Tipo_1;	
			InicializaFloat2D_MP(&MP_Corr2D_Tipo_1, NDX, NDY, 0);
			
Dist_MP MP_Histo_1;
float TamParticion1=1.0;
	InicializaDist_MP(&MP_Histo_1, TamParticion1, -1500,1000);
	

char contenedor[150];
		sprintf(contenedor, "../OBSERVACIONES");
//	CreaContenedor(contenedor);
	
/////////////////////////////////////Termina Prepara CONTENEDOR para escribir DATOS

	estado e;
	
	CargaDATOS("livetrees.txt",&e,NDX,NDY,TamParticion);

	ResetFloat2D_MP(&MP_Corr2D_Tipo_1);						
	CFFT_Tipos_MP(&e, 1, &MP_Corr2D_Tipo_1, 1, 1);
	ResetFloat1D_MP(&MP_Correlacion_1);
	CompactaCorrelacion(&MP_Corr2D_Tipo_1, &MP_Correlacion_1);		
	GuardaCorrelacion_MP(contenedor, "Auto", &MP_Correlacion_1);
	//GuardaCorrXY(&MP_Corr2D_1,"despues");

return;
}
