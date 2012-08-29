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
#include "ControlDinamico.h"

main(){	
int clave=124;
int *servicio = activa_escucha(clave);
//int *servicio;
//*servicio=0;
///////////////////////////Inicializa parametros de la simulacion
int NDX=500;
float fracON = 1.0;
int NDY=NDX;
int T_max = 10000;
int NoEnsambles=4;

float Birth1= 0.0;
float Coagulation1= 2.0; //Brown usa: 0.00002; 
float CoaIntra= 0.0000; //Modelo J-C 0.0008
float Dead1= 0.0;
int RadioBirth1= 10;
int RadioCoa1= 10;
int RadioCoaIntra1= 0;  //Modelo Heteromyopia 20

int CantidadEspecies=1;

float Birth2= 2.0;
float Coagulation2= 0; //Brown usa: 0.00002; 
float Dead2= 1.0;
int RadioBirth2= 10;
int RadioCoa2= 5;



//SetBirth(3.0,0);
//SetCoagulation(1.0,0);
//SetDead(1.0,0);
//SetCoagulationIntra(0.0,0);
//SetRadioBirth(RadioBirth1,NoEspecie);
//SetRadioCoa(RadioCoa1,NoEspecie);
//EscalaTiempoMetabolico(0);

int CantEspecies;
	for(CantEspecies=0;CantEspecies<=CantidadEspecies;CantEspecies++)
	{
			SetSpecie2(CantEspecies, Birth1, Coagulation1, CoaIntra, Dead1, RadioBirth1, RadioCoa1, RadioCoaIntra1);
	}

//	SetSpecie(1, Birth1, Coagulation1, Dead1, RadioBirth1, RadioCoa1);
//	SetSpecie(2, Birth2, Coagulation2, Dead2, RadioBirth2, RadioCoa2);

omp_set_num_threads(4);
////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:


Float2D_MP MP_RhoVsT_1;		
	int MaximoTamano = 500;
		InicializaFloat2D_MP(&MP_RhoVsT_1, T_max, MaximoTamano + 1, NoEnsambles);
		
		
Dist_MP MP_RhoDist_1;
float TamParticion=0.0001;
	InicializaDist_MP(&MP_RhoDist_1, TamParticion);
	
Float1D_MP TamDist_1;
		InicializaFloat1D_MP(&TamDist_1, T_max*2);
		
Float1D_MP TamDistRad_1;
		InicializaFloat1D_MP(&TamDistRad_1, T_max*2);
	
Float1D_MP MP_Correlacion_1;
	InicializaFloat1D_MP(&MP_Correlacion_1, NDX);

char contenedor[150];
	sprintf(contenedor,"SupCrecimiento0(B,D,C,RB,RC)@(%1.3f,%1.3f,%1.4f,%d,%d)_(NDX,Tmax)@(%d,%d)",Birth1,Dead1,Coagulation1,RadioBirth1,RadioCoa1,NDX,T_max);
	CreaContenedor(contenedor);
	
//char contenedorLec[150];	
	//sprintf(contenedorLec,"BrwRemNich0MP_(B,D,C,RB,RC)@(%1.3f,%1.3f,%1.3f,%d,%d)_(NDX,Tmax)@(%d,%d)",Birth1,Dead1,Coagulation1,RadioBirth1,RadioCoa1,NDX,T_max);
/////////////////////////////////////Termina Prepara CONTENEDOR para escribir DATOS



///////////////////////////////////// Estado INICIAL:

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
	
	int NoEspecie;
	for(NoEspecie=1;NoEspecie<=CantidadEspecies;NoEspecie++)
	{
		for(Par=0;Par<MaxPar;Par++)
		{
		GeneraEstadoAleatorioTamano(&e[Par], 0.5, 1, 1);
		GeneraEstadoAleatorioTamano(&e[Par], 1.0, -1, -1);
		}
	}

	
/////////////////////////////////////Termina Estado INICIAL
//////////////////////////Prepara Contenedor en Memoria de cada proceso Para Mejorar rendimiento (optimizar el uso de cache de cada procesador)

	Float2D_MP MP_RhoVsT;	
			InicializaFloat2D_MP(&MP_RhoVsT, T_max, MaximoTamano + 2, MaxPar);
			
	Dist_MP MP_RhoDist;
		InicializaDist_MP(&MP_RhoDist, TamParticion);
		
		Float1D_MP TamDist;
		InicializaFloat1D_MP(&TamDist, T_max*2);
		
		Float2D_MP MP_Corr2D_Tamano;
		InicializaFloat2D_MP(&MP_Corr2D_Tamano, NDX, NDY, 0);
		
		Float1D_MP MP_Correlacion;
		InicializaFloat1D_MP(&MP_Correlacion, NDX);
		
/////////////////////////////////Termina prepara Contenedor MEMORIA de cada PROCESO

		for(Par=0;Par<MaxPar;Par++)
		{
			ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
			ActualizaRecursos_MP(&e[Par],&MP_RhoVsT);	
		}
		
//////////////////////////////Barrido Monte CARLO:
	char distT[50];
	char prefix[10];
	int i;
	for(i=1;i<=T_max;i++)
	{
		for(Par=0;Par<MaxPar;Par++)
		{
			BarrMCcRyCampTamano(&e[Par], 1.0);
			ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
			ActualizaRecursos_MP(&e[Par],&MP_RhoVsT);	
			//#pragma omp master
			//{
				//GuardaEstadoEn(contenedor,&e[0]);
			//}
		}
		
		if(*servicio==1)    //Servicio de control Dinamico
		{
			printf("entro a servicio\n");
			Servicio(i,contenedor,&MP_RhoVsT, &MP_RhoVsT_1, &MP_RhoDist, &MP_RhoDist_1, e, MaxPar);
			SalidaCD(&i,T_max);
		} 
	
			if((i-(i/100)*100)==1)    //Inicializa cada 500 pasos
			{
				
				//printf("T:%d MT: %f, MM:%f\n",i,e[0].Meta_T,e[0].Max_Metabolic);
					for(Par=0;Par<MaxPar;Par++)
					{
						ActualizaDistTamano_MP(&e[Par], &TamDist_1, 'A');
					}
						
					//ResetFloat2D_MP(&MP_Corr2D_Tamano);						
					//CFFT_Mark_MP(e, MaxPar, &MP_Corr2D_Tamano, 0, 0);
					//ResetFloat1D_MP(&MP_Correlacion);
					//CompactaCorrelacion(&MP_Corr2D_Tamano, &MP_Correlacion);
					//SumaFloat1D_MP(&MP_Correlacion,&MP_Correlacion_1);
					
					#pragma omp barrier
					#pragma omp single
					{
						sprintf(distT,"DistTamanos_%d",i);
						GuardaFloat1D_MP(contenedor,distT,&TamDist_1);
						//sprintf(distT,"DistTamRad_%d",i);
						//GuardaFloat1D_MP(contenedor,distT,&TamDistRad_1);
						//sprintf(prefix,"Mark");
						//GuardaCorrelacion_MP(contenedor, prefix , &MP_Correlacion_1);
					}
			
				init_JKISS();
				//printf("Se ha reinicializado JKISS (c/100pasos) paso:%d, threath:%d \n",i,id);
			}

	}
	
//////////////////////////////Termina Monte CARLO

	for(Par=0;Par<MaxPar;Par++)
	{
		ActualizaDistTamano_MP(&e[Par], &TamDist_1, 'A');
	}
	SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
	//SumaFloat1D_MP(&TamDist,&TamDist_1);
	
	ResetFloat2D_MP(&MP_Corr2D_Tamano);						
	CFFT_Mark_MP(e, MaxPar, &MP_Corr2D_Tamano, 10, 10);
	ResetFloat1D_MP(&MP_Correlacion);
	CompactaCorrelacion(&MP_Corr2D_Tamano, &MP_Correlacion);
	SumaFloat1D_MP(&MP_Correlacion,&MP_Correlacion_1);
	sprintf(prefix,"T_%d_Mark",i);
	GuardaCorrelacion_MP(contenedor, prefix , &MP_Correlacion_1);
	
	
}	/////TERMINA PARALLEL

puts("Guardando informacion...\n");
GuardaRhoVsT_MP(contenedor,&MP_RhoVsT_1,NULL);
GuardaFloat1D_MP(contenedor,"DistTamanos",&TamDist_1);
//GuardaFloat1D_MP(contenedor,"DistTamRad",&TamDist_1);


cierra_escucha();

return;
}
