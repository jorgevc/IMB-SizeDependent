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

omp_set_num_threads(4);
////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:


//Float2D_MP MP_RhoVsT_1;
//int NoEspecies=CantidadEspecies;		
		//InicializaFloat2D_MP(&MP_RhoVsT_1, T_max, NoEspecies, NoEnsambles);
		
//Dist_MP MP_RhoDist_1;
//float TamParticion=0.0001;
//float Xini=0.0;
//float Xfin=0.12;
	//InicializaDist_MP(&MP_RhoDist_1, TamParticion,Xini,Xfin);
	
//Float1D_MP MP_Correlacion_1;
	//InicializaFloat1D_MP(&MP_Correlacion_1, NDX);
	
//Float2D_MP MP_Corr2D_1;	
		//	InicializaFloat2D_MP(&MP_Corr2D_1, NDX, NDY, 0);
			
//Float2D_MP MP_Corr2D_Tipo_1;	
			//InicializaFloat2D_MP(&MP_Corr2D_Tipo_1, NDX, NDY, 0);
			
Dist_MP MP_Histo_1;
float TamParticion1=1.0;
	InicializaDist_MP(&MP_Histo_1, TamParticion1, -1500,1000);
	

char contenedor[150];
	//sprintf(contenedor,"BrwRemJC50-1_(B,D,C,RB,RC)@(%1.3f,%1.3f,%1.4f:%1.4f,%d,%d:%d)_(NDX,Tmax)@(%d,%d)",Birth1,Dead1,Coagulation1,CoaIntra,RadioBirth1,RadioCoa1,RadioCoaIntra1,NDX,T_max);
//	sprintf(contenedor, "PruebasCorrelacion");
//	CreaContenedor(contenedor);
	
char contenedorLec[150];	
	sprintf(contenedorLec,"BrwRemHMRevised-1_(B,D,C,RB,RC)@(2.000,1.000,0.0000:0.0000,10,5:20)_(NDX,Tmax)@(100,100)");
	sprintf(contenedor,contenedorLec);
	
	
char tiempo[10];
	sprintf(tiempo,"T_087");	
/////////////////////////////////////Termina Prepara CONTENEDOR para escribir DATOS



///////////////////////////////////// INICIA PARALLEL

#pragma omp parallel			/////// Estado INICIAL:
{
	
	int num_threads = omp_get_num_threads();
	int id = omp_get_thread_num();
	int MaxPar = NoEnsambles/num_threads;
	#pragma omp master
	{
			MaxPar+= NoEnsambles - MaxPar * num_threads;	 
	}
	
	estado e[MaxPar];
	
	
	MaxPar=CargaEstado_MP(contenedorLec,tiempo,e,NDX,NDY,id,MaxPar);
	
	init_JKISS(); //Inicializa la semilla de cada proceso.
/////////////////////////////////////Termina Estado INICIAL
//////////////////////////Prepara Contenedor en Memoria de cada proceso Para Mejorar rendimiento (optimizar el uso de cache de cada nucleo)

	//Float2D_MP MP_RhoVsT;	
	//		InicializaFloat2D_MP(&MP_RhoVsT, T_max, NoEspecies, MaxPar);
			
//	Dist_MP MP_RhoDist;
	//	InicializaDist_MP(&MP_RhoDist, TamParticion,Xini,Xfin);

	Float2D_MP MP_Corr2D_Tipo;
		InicializaFloat2D_MP(&MP_Corr2D_Tipo, NDX, NDY, 0);
		
		Float1D_MP MP_Correlacion;
	InicializaFloat1D_MP(&MP_Correlacion, NDX);
	
		
/////////////////////////////////Termina prepara Contenedor MEMORIA de cada PROCESO
/*
		for(Par=0;Par<MaxPar;Par++)
		{
			ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);	
		}
		
//////////////////////////////Barrido Monte CARLO:
	int i;
	for(i=1;i<=T_max;i++)
	{
		for(Par=0;Par<MaxPar;Par++)
		{
			BarrMCcRyCamp(&e[Par]);
			ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
		}
		
		if(*servicio==1)    //Servicio de control Dinamico
		{
			printf("entro a servicio\n");
			Servicio(i,contenedor,&MP_RhoVsT, &MP_RhoVsT_1, &MP_RhoDist, &MP_RhoDist_1, e, MaxPar);
			SalidaCD(&i,T_max);
		} 
	
			if((i-(i/500)*500)==499)    //Inicializa cada 500 pasos
			{
				init_JKISS();
				printf("Se ha reinicializado JKISS (c/100pasos) paso:%d, threath:%d \n",i,id);
			}

	}
	
//////////////////////////////Termina Monte CARLO
*/
int Par;
	//for(Par=0;Par<MaxPar;Par++)
	//{
	//	ActualizaRhoVsT_MP(&e[Par],NULL,&MP_RhoDist);
	//	GuardaEstadoEn_MP(contenedor,&e[Par],id,Par);   //Guarda el ultimo estado de la corrida	
	//	ActualizaCorrelacion_MP(&e[Par],&MP_Correlacion);
	//	ActualizaCorrelacionTipo_MP(&e[Par], &MP_CorrelacionTipo, 1, 2);
	//}	
	

	//SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
	//SumaDist_MP(&MP_RhoDist,&MP_RhoDist_1);
	//SumaFloat1D_MP(&MP_Correlacion,&MP_Correlacion_1);
	//SumaFloat1D_MP(&MP_CorrelacionTipo,&MP_CorrelacionTipo_1);


#pragma omp barrier
	
	
	//CFFT_MP(e, MaxPar, &MP_Corr2D_1);
	int tipo1 = 1;
	int tipo2, NoPart;
	float Valor;
	char prefix[20];
int ENS = MaxPar;

	for(ENS=0;ENS<MaxPar;ENS+=10)
	{
		for(tipo1=1;tipo1<50;tipo1++)
		{
			for(tipo2 = tipo1 + 1; tipo2<=50;tipo2++)
			{
					ResetFloat2D_MP(&MP_Corr2D_Tipo);						
					CFFT_Tipos_MP(&e[ENS], 10, &MP_Corr2D_Tipo, tipo1, tipo2);
					if(MP_Corr2D_Tipo.NoEnsambles > 2)
					{
						ResetFloat1D_MP(&MP_Correlacion);
						CompactaCorrelacion(&MP_Corr2D_Tipo, &MP_Correlacion);
						//sprintf(prefix,"Id:%d_(%d,%d)",id,tipo1,tipo2);
						//GuardaCorrelacion_MP(contenedor, prefix , &MP_Correlacion);
						//SumaFloat1D_MP(&MP_Correlacion,&MP_Correlacion_1);
						//#pragma omp barrier
						//#pragma omp single
						//{
							//printf("Valor Corr %f\n",MP_Correlacion.array[50]);
							
							Valor=Integra(&MP_Correlacion,1,50);
							NoPart=(int)((Valor- MP_Histo_1.xIni)/MP_Histo_1.TamParticion);
							if(NoPart>=0 && NoPart<=MP_Histo_1.i_max)
							{
								MP_Histo_1.array[NoPart]++;
								MP_Histo_1.NoEnsambles++;
							}else{
								printf("No cupo %f de (%d,%d)\n",Valor,tipo1,tipo2);
							}
							//ResetFloat1D_MP(&MP_Correlacion_1);
						//}
					}
					//printf("Tipo %d vs %d echo\n",tipo1,tipo2, ENS);
				}
			}
			
			printf("ENS %d echo\n", ENS);
			//#pragma omp barrier
			//#pragma omp single
			//{
				//sprintf(nombre,"NonSelf_Area_Ensambles:%d",ENS);
				//GuardaDist_MP(contenedor,nombre,&MP_Histo_1);
				//printf("%s\n",nombre);
				//ResetDist_MP(&MP_Histo_1);
			//}
		}
			
		
	//GuardaCorrelacion_MP(contenedor,&MP_Correlacion);
	
}	/////TERMINA PARALLEL

/*
//GuardaCorrXY(&MP_Corr2D_1,"antes");
printf("Compactando Correlacion ...\n");
CompactaCorrelacion(&MP_Corr2D_1, &MP_Correlacion_1);
//GuardaCorrXY(&MP_Corr2D_1,"despues");
printf("Guadando Correlacion...\n");
GuardaCorrelacion_MP(contenedor,&MP_Correlacion_1);
ResetFloat1D_MP(&MP_Correlacion_1);
printf("Compactando Correlacion Tipos...\n");
CompactaCorrelacion(&MP_Corr2D_Tipo_1, &MP_Correlacion_1);
printf("Guadando Correlacion Tipos...\n");
GuardaCorrelacionTipo_MP(contenedor,&MP_Correlacion_1);


GuardaTiposEn_MP(contenedor,&MP_RhoVsT_1,T_max);
GuardaCorrelacion_MP(contenedor,&MP_Correlacion_1);
GuardaCorrelacionTipo_MP(contenedor,&MP_CorrelacionTipo_1); 
*/
puts("Guardando informacion...\n");
//GuardaRhoVsT_MP(contenedor,NULL,&MP_RhoDist_1);
char nombre[50];
sprintf(nombre,"xPOD_%s",tiempo);
GuardaDist_MP(contenedor,nombre,&MP_Histo_1);

//printf("Guadando Correlacion Tipos...\n");
//GuardaCorrelacion_MP(contenedor, "Auto", &MP_Correlacion_1);
 //printf("Guardando HistoGram\n");
//GuardaDist_MP(contenedor,"NonSelf_Area",&MP_Histo_1);
return;
}
