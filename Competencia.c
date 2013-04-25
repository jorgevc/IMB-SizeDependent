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
#include <math.h>


main(){	

///////////////////////////Inicializa parametros de la simulacion
runDescriptor run;
run.X=500;
run.Y=run.X;
run.grid_units=3.0;
run.size_units=9.0;

run.T_max=(run.grid_units*run.grid_units)*3000;
run.NoEnsambles=20;

float const Area_units=run.grid_units*run.grid_units;
float const Length_units=run.grid_units;

run.Model.coagulation_exp=1.0;  //no cambiar aqui
run.Model.coagulation_factor=(Area_units/pow(run.size_units,run.Model.coagulation_exp))*10.0; 
run.Model.coagulation_radio_exp=0.5; //no cambiar aqui
run.Model.coagulation_radio_factor=(Length_units/pow(run.size_units,run.Model.coagulation_radio_exp))*1.0;
run.Model.metabolic_exp=2.0;
run.Model.metabolic_factor=(Area_units/pow(run.size_units,run.Model.metabolic_exp))*0.01; 
run.Model.health_factor=0.1; //usandose lineal proporcional al tamano (adimensional) fraccion de biomasa que puede "danarse" antes de enfermar. 
run.Model.min_health=0;
run.resource_rate=1.0;
//(Area_units>=size_units)
run.Model.growth_constant=(Area_units/run.size_units)*100.0; // (int) needed resources per unit size increse.

run.Model.birth_rate=0.0;
run.Model.dead_rate=0.0;
run.Model.intra_coagulation=0.0;

int const write_interval=(run.grid_units*run.grid_units)*100;
int const separation=run.grid_units*80;
////
int T_max = run.T_max;
int NoEnsambles=run.NoEnsambles;


model modelo;
modelo = run.Model;

Individual indv;
indv.species=1;
indv.size=10;
indv.radio=1;
indv.metabolism=0;
indv.health=0;

int NDX=(run.X*run.grid_units)+1;
int NDY=(run.Y*run.grid_units)+1;
float delta_s=1.0/(run.size_units);

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
	sprintf(contenedor,"DATOS_TAM/25_April/4");
	CreaContenedor(contenedor,run);
	
Float1D_MP meanDensity;
	InicializaFloat1D_MP(&meanDensity, T_max+10);
	
Float1D_MP meanSize;
	InicializaFloat1D_MP(&meanSize, T_max+10);

FILE *file;

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
			e[Par].units=run.grid_units;
			e[Par].size_units=run.size_units;
		}

		init_JKISS(); //Inicializa la semilla de cada proceso.
		
			
		int i,j;
			for(Par=0;Par<MaxPar;Par++)
			{
			//GeneraEstadoAleatorioTamano(&e[Par], 0.5, 1, 1);
				for(i=separation;i<NDX;i+=separation)
				{
					for(j=separation;j<NDY;j+=separation)
					{
						InsertIndividualAt(&e[Par],i,j,indv,1);
					}
				}
			GeneraEstadoAleatorioTamano(&e[Par], run.resource_rate, -1, -1);
			setMaxMetabolic(&e[Par],&modelo);
			}

		
	/////////////////////////////////////Termina Estado INICIAL
	//////////////////////////Prepara Contenedor en Memoria de cada proceso Para Mejorar rendimiento (optimizar el uso de cache de cada procesador)

		Float2D_MP MP_RhoVsT;	
				InicializaFloat2D_MP(&MP_RhoVsT, T_max, MaximoTamano, MaxPar);

		Float1D_MP TamDist;
		InicializaFloat1D_MP(&TamDist, T_max+10);
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
		
		Rate_log rate[5];
		int MaximoTamanoIni = 100;
		InitRate_log(&rate[0],MaximoTamanoIni);
		InitRate_log(&rate[1],MaximoTamanoIni);
		InitRate_log(&rate[2],MaximoTamanoIni);
		InitRate_log(&rate[3],MaximoTamanoIni);
		InitRate_log(&rate[4],MaximoTamanoIni);
		FILE *taza;
		int r;
		
		//char prefix[10];
		//int i;
		for(i=1;i<=T_max;i++)
		{		
			for(Par=0;Par<MaxPar;Par++)
			{
				BarrMCcRyCampTamano(&e[Par], run.resource_rate, &modelo, rate);
				ActualizaRecursos_MP(&e[Par],&MP_RhoVsT);
				ActualizaDistTamano_MP(&e[Par], &TamDist, 'A');
				#pragma omp atomic
				meanDensity.array[e[Par].T]+=e[Par].ON;
			}
			
			#pragma omp barrier
			#pragma omp atomic
			meanSize.array[TamDist.T]+=TamDist.array[0];
			
				if((i-(i/write_interval)*write_interval)==1)    //Inicializa cada write_interval
				{
					//for(Par=0;Par<MaxPar;Par++)
					//{
					//	ActualizaDistTamano_MP(&e[Par], &TamDist_1, 'A');
					//}
					
					SumaFloat1D_MP(&TamDist,&TamDist_1);
					#pragma omp barrier
					#pragma omp single
					{
						sprintf(distT,"DT_%d",i);	
						GuardaFloat1D_MP(contenedor,distT,&TamDist_1);
						ResetFloat1D_MP(&TamDist_1);
					}
					#pragma omp master
					{
						GuardaEstadoEn(contenedor, &e[0]);
					}
					init_JKISS();
				}

		}
		
	//////////////////////////////Termina Monte CARLO
	#pragma omp barrier
	#pragma omp single
	{
		sprintf(distT,"%s/GrowthR",contenedor,i);
		taza=fopen(distT, "w");						
		for(r=1;r<=rate[0].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*rate[0].Growth[r]/((float)rate[0].NoEnsambles[r]));
		}
		fclose(taza);
		
		sprintf(distT,"%s/DeadR",contenedor,i);
		taza=fopen(distT, "w");
		for(r=1;r<=rate[1].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, rate[1].Growth[r]/((float)rate[1].NoEnsambles[r]));
		}
		fclose(taza);
		
		sprintf(distT,"%s/ResourceR",contenedor,i);
		taza=fopen(distT, "w");
		for(r=1;r<=rate[2].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s,delta_s*rate[2].Growth[r]/((float)rate[2].NoEnsambles[r]));
		}
		fclose(taza);
		
		sprintf(distT,"%s/MetabolicR",contenedor,i);
		taza=fopen(distT, "w");
		for(r=1;r<=rate[3].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*rate[3].Growth[r]/((float)rate[3].NoEnsambles[r]));
		}
		fclose(taza);
	
	}

		SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
		
		//Libera Memoria
		for(Par=0;Par<MaxPar;Par++)
		{
			LiberaMemoria(&e[Par]);
		}
		
		LiberaMemoriaFloat2D_MP(&MP_RhoVsT);
	}	/////TERMINA PARALLEL

GuardaRhoVsT_MP(contenedor,&MP_RhoVsT_1,NULL);
LiberaMemoriaFloat2D_MP(&MP_RhoVsT_1);
LiberaMemoriaFloat1D_MP(&TamDist_1);

char thinning[50];
sprintf(thinning,"%s/thinning",contenedor);
file=fopen(thinning, "w");
fputs("# T meanSize meanDensity\n",file);
int r;
for(r=1;r<=T_max;r++){
fprintf(file,"%d %f %f\n",r, delta_s*(meanSize.array[r]/NoEnsambles), meanDensity.array[r]/NoEnsambles);
}
fclose(file);

return;
}
