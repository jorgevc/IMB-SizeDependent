#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "libPP_5.0.h"
#include "EntSalArb_MP.h"
#include "GNA.h"
#include "ControlDinamico.h"

main(){	
int clave=323;
int *servicio = activa_escucha(clave);
//int *servicio;
//*servicio=0;
///////////////////////////Inicializa parametros de la simulacion
int NDX=100;
float fracON = 1.0;
int NDY=NDX;
int T_max = 1000;
int NoEnsambles=1000;

float Birth1= 2.0;
float Coagulation1= 0.00002; //Brown usa: 0.00002; 
float CoaIntra= 0.00002; //Modelo J-C 0.0008
float Dead1= 1.0;
int RadioBirth1= 10;
int RadioCoa1= 5;
int RadioCoaIntra1= 10;  //Modelo Heteromyopia 20

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

omp_set_num_threads(4);
////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:


Float2D_MP MP_RhoVsT_1;
int NoEspecies=CantidadEspecies;		
		InicializaFloat2D_MP(&MP_RhoVsT_1, T_max, NoEspecies, 0);
		
Dist_MP MP_RhoDist_1;
float TamParticion=0.0001;
float Xini=0.0;
float Xfin=0.2;
	InicializaDist_MP(&MP_RhoDist_1, TamParticion,Xini,Xfin);
	
Float1D_MP MP_Correlacion_1;
	InicializaFloat1D_MP(&MP_Correlacion_1, NDX);
	
Float2D_MP MP_Corr2D_1;
	InicializaFloat2D_MP(&MP_Corr2D_1, NDX, NDY, 0);
			
//Float2D_MP MP_Corr2D_Tipo_1;	
			//InicializaFloat2D_MP(&MP_Corr2D_Tipo_1, NDX, NDY, 0);

char contenedor[150];
	sprintf(contenedor,"BrwRemNeutralRevised-1_(B,D,C,RB,RC)@(%1.3f,%1.3f,%1.4f:%1.4f,%d,%d:%d)_(NDX,Tmax)@(%d,%d)",Birth1,Dead1,Coagulation1,CoaIntra,RadioBirth1,RadioCoa1,RadioCoaIntra1,NDX,T_max);
	//sprintf(contenedor,"BrwRemMPHM-2_(B,D,C,RB,RC)@(2.000,1.000,0.000,10,5)_(NDX,Tmax)@(1000,150)");
	CreaContenedor(contenedor);
	
char contenedorLec[150];	
	//sprintf(contenedorLec,"BrwRemNich0MP_(B,D,C,RB,RC)@(%1.3f,%1.3f,%1.3f,%d,%d)_(NDX,Tmax)@(%d,%d)",Birth1,Dead1,Coagulation1,RadioBirth1,RadioCoa1,NDX,T_max);
	sprintf(contenedorLec,contenedor);
/////////////////////////////////////Termina Prepara CONTENEDOR para escribir DATOS



///////////////////////////////////// INICIA PARALLEL

#pragma omp parallel			///////Estado INICIAL:
{
	init_JKISS(); //Inicializa la semilla de cada proceso.
	
	int num_threads = omp_get_num_threads();
	int id = omp_get_thread_num();
	int MaxPar = NoEnsambles/num_threads;
	#pragma omp master
	{
			MaxPar+= NoEnsambles - MaxPar * num_threads;	 
	}
	estado e[MaxPar];
	
	//MaxPar=CargaEstado_MP(contenedorLec,"T_001",e,NDX,NDY,id,MaxPar);
	
	int Par;
	for(Par=0;Par<MaxPar;Par++)
	{
		AlojaMemoria(&e[Par], NDX, NDY);
		ResetEstado(&e[Par]); 	
	}

	
	int NoEspecie;
	for(NoEspecie=1;NoEspecie<=CantidadEspecies;NoEspecie++)
	{
		for(Par=0;Par<MaxPar;Par++)
		{
		//InsertaIndividuosAleatorio(&e[Par],100,NoEspecie);
		GeneraEstadoAleatorio(&e[Par], 0.007, NoEspecie);
		}
	}

	
/////////////////////////////////////Termina Estado INICIAL
//////////////////////////Prepara Contenedor en Memoria de cada proceso Para Mejorar rendimiento (optimizar el uso de cache de cada procesador)

	Float2D_MP MP_RhoVsT;	
			InicializaFloat2D_MP(&MP_RhoVsT, T_max, NoEspecies, MaxPar);
			
	Dist_MP MP_RhoDist;
		InicializaDist_MP(&MP_RhoDist, TamParticion,Xini,Xfin);
		
		
/////////////////////////////////Termina prepara Contenedor MEMORIA de cada PROCESO

		//for(Par=0;Par<MaxPar;Par++)
		//{
			//ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);	
		//}
		
//////////////////////////////Barrido Monte CARLO:
	int i;
	for(i=0;i<T_max;i++)
	{
		for(Par=0;Par<MaxPar;Par++)
		{
			if((i-(i/10)*10)==0) // Guarda Estado cada 10 pasos
			{
			GuardaEstadoEn_MP(contenedor,&e[Par],id,Par);
			}
			BarrMCcRyCamp(&e[Par]);
			ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
		}
		
		if(*servicio==1)    //Servicio de control Dinamico
		{
			printf("entro a servicio\n");
			Servicio(i,contenedor,&MP_RhoVsT, &MP_RhoVsT_1, &MP_RhoDist, &MP_RhoDist_1, e, MaxPar);
			SalidaCD(&i,T_max);
		} 
	
			if((i-(i/500)*500)==499)    //Inicializa cada 1000 pasos
			{
				init_JKISS();
				printf("Se ha reinicializado JKISS (c/100pasos) paso:%d, threath:%d \n",i,id);
			}

	}
	
//////////////////////////////Termina Monte CARLO
	printf("Guardando ultimo estado...\n");
	for(Par=0;Par<MaxPar;Par++)
	{
		//ActualizaRhoVsT_MP(&e[Par],NULL,&MP_RhoDist);
		GuardaEstadoEn_MP(contenedor,&e[Par],id,Par);   //Guarda el ultimo estado de la corrida	
	}	
	

	SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
	//SumaDist_MP(&MP_RhoDist,&MP_RhoDist_1);
	//printf("Calculando correlacion...\n");
	//CFFT_MP(e, MaxPar, &MP_Corr2D_1);
	//CFFT_Tipos_MP(e, MaxPar, &MP_Corr2D_1, 5, 5);
		
}	////////////////////////////////////////////////////////////////////TERMINA PARALLEL

//puts("Guardando informacion...\n");
GuardaRhoVsT_MP(contenedor,&MP_RhoVsT_1,NULL);
//GuardaTiposEn_MP(contenedor,&MP_RhoVsT_1,T_max);
//printf("Compactando Correlacion ...\n");
//CompactaCorrelacion(&MP_Corr2D_1, &MP_Correlacion_1);
//printf("Guadando Correlacion...\n");
//GuardaCorrelacion_MP(contenedor, "Auto" , &MP_Correlacion_1);
//GuardaCorrelacion_MP(contenedor,&MP_Correlacion_1);
//ResetFloat1D_MP(&MP_Correlacion_1);
//printf("Compactando Correlacion tipo...\n");
//CompactaCorrelacion(&MP_Corr2D_Tipo_1, &MP_Correlacion_1);
//printf("Guadando Correlacion tipo...\n");
//GuardaCorrelacionTipo_MP(contenedor,&MP_Correlacion_1);


cierra_escucha();

return;
}
