/*
Copyright 2012 Jorge Velazquez
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "libPP_6.1.h"
#include "EntSalArb_MP_Comp.h"
#include "GNA.h"

main(){	

///////////////////////////Inicializa parametros de la simulacion
int NDX=500;
int NDY=NDX;
int T_max = 1500;
int NoEnsambles=4;

float Birth1= 0.0;
float Coagulation1= 2.0; //Brown usa: 0.00002; 
float CoaIntra= 0.0000; //Modelo J-C 0.0008
float Dead1= 0.0;
int RadioBirth1= 10;
int RadioCoa1= 10;
int RadioCoaIntra1= 0;  //Modelo Heteromyopia 20

int CantidadEspecies=1;

model modelo;

int CantEspecies;
	for(CantEspecies=0;CantEspecies<=CantidadEspecies;CantEspecies++)
	{
			SetSpecie2(CantEspecies, Birth1, Coagulation1, CoaIntra, Dead1, RadioBirth1, RadioCoa1, RadioCoaIntra1);
	}

omp_set_num_threads(4);
////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:



Float2D_MP MP_RhoVsT_1;		
	int MaximoTamano = 1;
		InicializaFloat2D_MP(&MP_RhoVsT_1, T_max, MaximoTamano, 0);
		
//Dist_MP MP_RhoDist_1;
//float TamParticion=0.0001;
//	InicializaDist_MP(&MP_RhoDist_1, TamParticion);
	
Float1D_MP TamDist_1;
		InicializaFloat1D_MP(&TamDist_1, T_max+10);
		
//Float1D_MP Experiment_1;
	//	InicializaFloat1D_MP(&Experiment_1, 41);

//Float1D_MP MP_Correlacion_1;
//	InicializaFloat1D_MP(&MP_Correlacion_1, NDX);

char contenedor[150];
	sprintf(contenedor,"DATOS_TAM/DI=0.5_CF=1.0_MF=0.5:50:1.5_RC=10_NX=500");
	CreaContenedor(contenedor);
	
	
//char contenedorLec[150];	
	//sprintf(contenedorLec,"BrwRemNich0MP_(B,D,C,RB,RC)@(%1.3f,%1.3f,%1.3f,%d,%d)_(NDX,Tmax)@(%d,%d)",Birth1,Dead1,Coagulation1,RadioBirth1,RadioCoa1,NDX,T_max);
/////////////////////////////////////Termina Prepara CONTENEDOR para escribir DATOS


///////////////////////////////////// Estado INICIAL:

modelo.CoaFact=1.0;
modelo.CoaExp=1.0;
modelo.MetFact=0.5;
modelo.MetExp=1.0;
modelo.ResurcesFact=1.0;


	#pragma omp parallel			///////INICIA PARALLEL
	{
		
		int num_threads = omp_get_num_threads();
		int id = omp_get_thread_num();
		int MaxPar = NoEnsambles/num_threads;
		#pragma omp master
		{
				MaxPar+= NoEnsambles - MaxPar * num_threads;	 
		}
		estado e[MaxPar];
		
		//MaxPar=CargaEstado_MP(contenedorLec,"T_8000",e,NDX,NDY,id,MaxPar);
		
		int Par;
		for(Par=0;Par<MaxPar;Par++)
		{
			AlojaMemoria(&e[Par], NDX, NDY);
			ResetEstado(&e[Par]); 	
		}

		init_JKISS(); //Inicializa la semilla de cada proceso.
		
			for(Par=0;Par<MaxPar;Par++)
			{
			GeneraEstadoAleatorioTamano(&e[Par], 0.5, 1, 1);
			GeneraEstadoAleatorioTamano(&e[Par], 1.0, -1, -1);
			}

		
	/////////////////////////////////////Termina Estado INICIAL
	//////////////////////////Prepara Contenedor en Memoria de cada proceso Para Mejorar rendimiento (optimizar el uso de cache de cada procesador)

		Float2D_MP MP_RhoVsT;	
				InicializaFloat2D_MP(&MP_RhoVsT, T_max, MaximoTamano, MaxPar);

		//Dist_MP MP_RhoDist;
		//	InicializaDist_MP(&MP_RhoDist, TamParticion);
			
		//	Float1D_MP TamDist;
		//	InicializaFloat1D_MP(&TamDist, T_max*2);
			
		//	Float2D_MP MP_Corr2D_Tamano;
		//	InicializaFloat2D_MP(&MP_Corr2D_Tamano, NDX, NDY, 0);
			
		//	Float1D_MP MP_Correlacion;
		//	InicializaFloat1D_MP(&MP_Correlacion, NDX);
			
	/////////////////////////////////Termina prepara Contenedor MEMORIA de cada PROCESO

			//for(Par=0;Par<MaxPar;Par++)
			//{
				//ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
				//ActualizaRecursos_MP(&e[Par],&MP_RhoVsT);	
			//}
			
	//////////////////////////////Barrido Monte CARLO:
		char distT[50];
		//char prefix[10];
		int i;
		for(i=1;i<=T_max;i++)
		{
			for(Par=0;Par<MaxPar;Par++)
			{
				BarrMCcRyCampTamano(&e[Par], 1.0, &modelo);
				//ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
				ActualizaRecursos_MP(&e[Par],&MP_RhoVsT);	
			}
							
				if((i-(i/50)*50)==1)    //Inicializa cada 100 pasos
				{
					
					for(Par=0;Par<MaxPar;Par++)
					{
						ActualizaDistTamano_MP(&e[Par], &TamDist_1, 'A');
					}
			
					#pragma omp barrier
					#pragma omp single
					{
						sprintf(distT,"DistTamanos_%d",i);
						GuardaFloat1D_MP(contenedor,distT,&TamDist_1);
					}
							
				
					init_JKISS();
				}

		}
		
	//////////////////////////////Termina Monte CARLO

		SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
		
		//Libera Memoria
		for(Par=0;Par<MaxPar;Par++)
		{
			LiberaMemoria(&e[Par]);
		}
		
		LiberaMemoriaFloat2D_MP(MP_RhoVsT);
	}	/////TERMINA PARALLEL

GuardaRhoVsT_MP(contenedor,&MP_RhoVsT_1,NULL);
LiberaMemoriaFloat2D_MP(MP_RhoVsT_1);
LiberaMemoriaFloat1D_MP(TamDist_1);
return;
}
