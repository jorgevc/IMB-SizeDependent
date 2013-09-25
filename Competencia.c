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
#ifdef VIRTUAL_GRID
printf("VIRUAL_GRID=TRUE\n");
#endif
#ifdef SOI
printf("SOI=TRUE\n");
#endif

///////////////////////////Inicializa parametros de la simulacion
runDescriptor run;
run.X=200;  	//unidades "fisicas"
run.Y=run.X;
run.grid_units=1.0; //factor de conversion de unidades "fisicas" a lado de celda
run.size_units=1.0; //numero de unidades en computo que hacen una unidad "fisica" de tamano (divisor de conversion)

//run.T_max=(run.grid_units*run.grid_units)*802;
run.T_max=(run.grid_units*run.grid_units)*20000;
//run.T_max=(run.grid_units*run.grid_units)*20000;
run.NoEnsambles=20;

//in this case the Area_units and Length_units are units of the imaginary resources grid.
run.Model.ResourcesScale = 50; //conversion factor to the imaginary resources grid or SOI resolution
run.Model.competitionAsymetry = 1.0;
double const Area_units=(run.grid_units*run.Model.ResourcesScale)*(run.grid_units*run.Model.ResourcesScale);
double const Length_units=run.grid_units*run.Model.ResourcesScale;

double coagulation_units = 1.0;
#ifdef SOI
coagulation_units = 3000.0;
run.Model.coagulation_factor = coagulation_units*1.0;
#else
run.Model.coagulation_factor = coagulation_units*1.0*3.1416;
#endif
//run.Model.coagulation_factor=((Area_units/pow(run.size_units,run.Model.coagulation_exp))*1.0);

run.Model.coagulation_radio_exp=0.375; //cambiar tambien en model.c  :0.375
run.Model.coagulation_exp=2*run.Model.coagulation_radio_exp;  // cambiar tambien en model.c
run.Model.coagulation_radio_factor=((Length_units)/pow(run.size_units,run.Model.coagulation_radio_exp))*1.0;
run.Model.metabolic_exp=1.0; //cambiar tambien en model.c
run.Model.metabolic_factor=coagulation_units*(Area_units/pow(run.size_units,run.Model.metabolic_exp))*0.75; 
// 0.2 anterior 1.25
run.Model.health_factor=2.0; //usandose lineal proporcional al tamano (adimensional) fraccion de biomasa que puede "danarse" antes de enfermar. 
// con 2.0 salen buenos resultados para U shape de dead rate.
run.Model.growth_constant=(coagulation_units*(Area_units/run.size_units)*0.05); // (int)>0 needed resources per unit size increse.
//run.Model.growth_constant=5;
run.Model.resource_rate=1.0; // morir por vejez
//run.Model.resource_rate=0.95;
run.initialMeanDistance=run.grid_units*12;
run.initialMinSeparation=10;

#ifdef HEALTH_TRACK	
run.Model.min_health=0;
#else
run.Model.min_health=0;
#endif


//run.Model.birth_rate=0.5;
run.Model.birth_rate=0.5;
run.Model.RadioBirth=10;
run.Model.dead_rate=0.0;
run.Model.intra_coagulation=0.0;

int const write_interval=(run.grid_units*run.grid_units)*200;
//int const separation=run.grid_units*3;
float const separation=(float)run.initialMeanDistance;
////
int T_max = run.T_max;
int NoEnsambles=run.NoEnsambles;


model modelo;
modelo = run.Model;


Individual indv;
indv.species=1;
indv.size=run.size_units*2;
indv.radio=R(indv, (&modelo));
indv.metabolism=0;
indv.health=0; 


int NDX=(run.X*run.grid_units)+1;
int NDY=(run.Y*run.grid_units)+1;
double delta_s=1.0/(run.size_units);

omp_set_num_threads(4);

////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:



//Float2D_MP MP_RhoVsT_1;		
//	int MaximoTamano = 1;
//		InicializaFloat2D_MP(&MP_RhoVsT_1, T_max+10, MaximoTamano, 0);
		
//Dist_MP MP_RhoDist_1;
//float TamParticion=0.0001;
//	InicializaDist_MP(&MP_RhoDist_1, TamParticion);
	
Float1D_MP TamDist_1;
		InicializaFloat1D_MP(&TamDist_1, T_max+10);
		
//Float1D_MP Experiment_1;
	//	InicializaFloat1D_MP(&Experiment_1, 41);

Float1D_MP MP_CorrelacionG;
	InicializaFloat1D_MP(&MP_CorrelacionG, NDX);
	ResetFloat1D_MP(&MP_CorrelacionG);
	

char contenedor[150];
	sprintf(contenedor,"DATOS_TAM/25_Sep/SteadyStateMaturity60Dead2Asym1Rand");
	CreaContenedor(contenedor,run);
	
Float1D_MP meanDensity;
	InicializaFloat1D_MP(&meanDensity, T_max+10);
	
Float1D_MP meanSize;
	InicializaFloat1D_MP(&meanSize, T_max+10);
	
float meanResourceG[T_max+10];
	
	
float time_map[T_max+10];
time_map[0]=0.0;

Rate_log Grate[6];
InitRate_log(&Grate[5],10);

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
			GeneraEstadoAleatorioTamano(&e[Par], 1.0/separation , indv);		
		//	FilterMinDistance(&e[Par], separation);
		FilterMinDistance(&e[Par], run.initialMinSeparation);
			
		//		for(i=separation;i<NDX;i+=separation)
			//	{
				//	for(j=separation;j<NDY;j+=separation)
					//{
						//InsertIndividualAt(&e[Par],i,j,indv,1);
					//}
				//}
			
			setMaxMetabolic(&e[Par],&modelo);
			
			}

		
	/////////////////////////////////////Termina Estado INICIAL
	//////////////////////////Prepara Contenedor en Memoria de cada proceso Para Mejorar rendimiento (optimizar el uso de cache de cada procesador)

//		Float2D_MP MP_RhoVsT;	
//				InicializaFloat2D_MP(&MP_RhoVsT, T_max, MaximoTamano, MaxPar);

		Float1D_MP TamDist;
		InicializaFloat1D_MP(&TamDist, T_max+100);
		//Dist_MP MP_RhoDist;
		//	InicializaDist_MP(&MP_RhoDist, TamParticion);
			
		//	Float1D_MP TamDist;
		//	InicializaFloat1D_MP(&TamDist, T_max*2);
			
			Float2D_MP MP_Corr2D;
			InicializaFloat2D_MP(&MP_Corr2D, NDX, NDY, 0);
			
			Float1D_MP MP_Correlacion;
			InicializaFloat1D_MP(&MP_Correlacion, NDX);
			
	/////////////////////////////////Termina prepara Contenedor MEMORIA de cada PROCESO

			//for(Par=0;Par<MaxPar;Par++)
			//{
				//ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
				//ActualizaRecursos_MP(&e[Par],&MP_RhoVsT);	
			//}
			
	//////////////////////////////Barrido Monte CARLO:
		char distT[50];
		char corrName[50];
		Rate_log rate[6];
		int MaximoTamanoIni = 100;
		InitRate_log(&rate[0],MaximoTamanoIni);
		InitRate_log(&rate[1],MaximoTamanoIni);
		InitRate_log(&rate[2],MaximoTamanoIni);
		InitRate_log(&rate[3],MaximoTamanoIni);
		InitRate_log(&rate[4],MaximoTamanoIni);
		InitRate_log(&rate[5],rate[2].i_max+100);
		
		Grupo originGroup,targetGroup;
		originGroup.TIPO=0; //all species
		originGroup.NEG=0; //lugares ocupados (contrario a lugares vacios)
		targetGroup.TIPO=0;
		targetGroup.NEG=0;
		CorrDescriptor corrEspecification;
		corrEspecification.MeanSquare=0;
		corrEspecification.NoEnsambles=MaxPar;
		corrEspecification.NoMuestras=1;
		corrEspecification.Muestra=0;
		Individual indvTmp,indvTmp2;
		indv.species=1;
		int s1,s2,r2;
		float Interaction,s1val,s2val,f, Area,Iss;
		FILE *fileKappa;
		FILE *fileCorr;
		
		for(i=1;i<=T_max;i++)
		{			
			for(Par=0;Par<MaxPar;Par++)
			{
				BarrMCcRyCampTamano(&e[Par], run.Model.resource_rate, &modelo, rate);			
				ActualizaDistTamano_MP(&e[Par], &TamDist, 'A');				
				#pragma omp atomic
				meanDensity.array[e[Par].T]+=e[Par].ON;
			}
			
			#pragma omp barrier
			#pragma omp atomic
			meanSize.array[TamDist.T]+=TamDist.array[0];
			#pragma omp master
			{
				time_map[e[0].T]=e[0].Meta_T;
			}
			
			
				if((i-(i/write_interval)*write_interval)==1)    //Inicializa cada write_interval
				{
					//for(Par=0;Par<MaxPar;Par++)
					//{
					//	ActualizaDistTamano_MP(&e[Par], &TamDist_1, 'A');
					//}
					
					///// meanK
					#pragma omp single
					{
						FreeRate_log(&Grate[5]);
						InitRate_log(&Grate[5],rate[5].i_max+200);
					}
					SumRate_log(&rate[5], &Grate[5]);
					FreeRate_log(&rate[5]);
					InitRate_log(&rate[5],rate[2].i_max);
					#pragma omp barrier
					#pragma omp master
					{			
						sprintf(distT,"%s/MeanK",contenedor);
								file=fopen(distT, "a");						
								for(j=1;j<=Grate[5].i_max;j++){
									if(Grate[5].NoEnsambles[j]>9)
									{
									fprintf(file,"%f %d %f  \n",((float)j)*delta_s, e[0].T , delta_s*Grate[5].Growth[j]/(coagulation_units*(float)Grate[5].NoEnsambles[j]));
									}
								}
								fclose(file);	
					}				
					//	
					//distribution
					#pragma omp single
					{
						ResetFloat1D_MP(&TamDist_1);
					}
					SumaFloat1D_MP(&TamDist,&TamDist_1);
					#pragma omp barrier
					#pragma omp single
					{
						sprintf(distT,"DT_%d",i);	
						GuardaFloat1D_MP(contenedor,distT,&TamDist_1);			
					}
					//
					//Analitical Mean Resorce Intake
					//#pragma omp master
					//{
					//sprintf(distT,"%s/kappa_Iss",contenedor);
					//fileKappa=fopen(distT, "a");
					//fputs("r1  r2  kappa(Integral)  T\n",fileKappa);
					
					//sprintf(corrName,"%s/corr", contenedor);
					//fileCorr=fopen(corrName, "a");
					//fputs("r1  r2  r  g  T\n",fileCorr);
					//}
					
					//for(s1=1;s1<=Grate[5].i_max;s1++)
					//{
						//if(TamDist_1.array[s1] > 0.0 )
						//{
							//indvTmp.size=run.size_units*s1;
							//indvTmp.radio=R(indvTmp, (&modelo));
							//Interaction=0.0;
							//for(s2=1;s2<=Grate[5].i_max;s2++)
							//{
								//ResetFloat2D_MP(&MP_Corr2D);
								//ResetFloat1D_MP(&MP_Correlacion);
								//originGroup.s=s1;
								//targetGroup.s=s2;
								//CFFT_Univ_MP(e, &corrEspecification, &MP_Corr2D, &originGroup, &targetGroup);
								//if(MP_Corr2D.NoEnsambles > 2)
								//{
									//CompactaCorrelacion(&MP_Corr2D, &MP_Correlacion);
									//SumaFloat1D_MP(&MP_Correlacion,&MP_CorrelacionG);
								//}
								//#pragma omp barrier
								//#pragma omp master
								//{
									//s1val=pow((double)s1 , run.Model.competitionAsymetry);
									//s2val=pow((double)s2 , run.Model.competitionAsymetry);
									//f = s1val/(s1val + s2val);	
									//indvTmp2.size=run.size_units*s2;
									//indvTmp2.radio=R(indvTmp2, (&modelo));
									//Iss=IntegraAC(&MP_CorrelacionG, indvTmp.radio,indvTmp2.radio,Length_units,i);
									//Interaction+=(TamDist_1.array[s2]*meanDensity.array[i]*(1-f)*Iss)/((float)(TamDist_1.NoEnsambles*NoEnsambles*NDX*NDY));
									//if(s1==s2)
									//{
										//for(r2=1;r2<=MP_CorrelacionG.i_max;r2++)
										//{
										//fprintf(fileCorr,"%d  %d  %d  %f  %d\n",indvTmp.radio,indvTmp2.radio, ((int)Length_units)*r2, MP_CorrelacionG.array[r2], i);
										//}
									//}		
										//if(Iss>0.0)
										//{
											//fprintf(fileKappa,"%d  %d  %f  %d\n", indvTmp.radio, indvTmp2.radio, Iss, i);
										//}
										
									//ResetFloat1D_MP(&MP_CorrelacionG);		
								//}
							//}
							//#pragma omp master
							//{
							//indvTmp.size=run.size_units*s1;
							//indvTmp.radio=R(indvTmp, (&modelo));
							//Area=3.1416*indvTmp.radio*indvTmp.radio;
							//meanResourceG[s1]=run.Model.resource_rate*(Area - Interaction);
							//}
						//}else{
							//meanResourceG[s1]=0.0;
						//}
					//}			
					//#pragma omp master
					//{
						//fclose(fileKappa);	
						//fclose(fileCorr);
							//sprintf(distT,"%s/analiticResource",contenedor);
								//file=fopen(distT, "a");
								//fputs("s  meanResource  T",file);
								//for(s1=1;s1<=Grate[5].i_max;s1++){
									//if(meanResourceG[s1]!=0.0)
									//{
										//fprintf(file,"%d  %f  %d\n", s1, meanResourceG[s1], i);
									//}
								//}
								//fclose(file);						
					//}
					//
					// estate
					#pragma omp master
					{
						GuardaEstadoEn(contenedor, &e[0]);
					}
					//
					init_JKISS();
				}

		}
		
	//////////////////////////////Termina Monte CARLO
	LiberaMemoriaFloat2D_MP(&MP_Corr2D);
	LiberaMemoriaFloat1D_MP(&MP_Correlacion);

	#pragma omp single
	{
		InitRate_log(&Grate[0],rate[0].i_max);
		InitRate_log(&Grate[1],rate[1].i_max);
		InitRate_log(&Grate[2],rate[2].i_max);
		InitRate_log(&Grate[3],rate[3].i_max);
		InitRate_log(&Grate[4],rate[4].i_max);
	}

	SumRate_log(&rate[0], &Grate[0]);
	SumRate_log(&rate[1], &Grate[1]);
	SumRate_log(&rate[2], &Grate[2]);
	SumRate_log(&rate[3], &Grate[3]);
	SumRate_log(&rate[4], &Grate[4]);

	FreeRate_log(&rate[0]);
	FreeRate_log(&rate[1]);
	FreeRate_log(&rate[2]);
	FreeRate_log(&rate[3]);
	FreeRate_log(&rate[4]);
	FreeRate_log(&rate[5]);
	
	//	#pragma omp master
	//	{
	//	SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
	//	}
		
		//Libera Memoria
		for(Par=0;Par<MaxPar;Par++)
		{
			LiberaMemoria(&e[Par]);
		}
		
	//	LiberaMemoriaFloat2D_MP(&MP_RhoVsT);
	}	/////TERMINA PARALLEL

// Write Rate
FILE *taza;
int r;
char distT[50];
		sprintf(distT,"%s/GrowthR",contenedor);
		taza=fopen(distT, "w");						
		for(r=1;r<=Grate[0].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*Grate[0].Growth[r]/((double)Grate[0].NoEnsambles[r]));
		}
		fclose(taza);
		
		sprintf(distT,"%s/DeadR",contenedor);
		taza=fopen(distT, "w");
		for(r=1;r<=Grate[1].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, Grate[1].Growth[r]/((float)Grate[1].NoEnsambles[r]));
		}
		fclose(taza);
		
		sprintf(distT,"%s/ResourceR",contenedor);
		taza=fopen(distT, "w");
		for(r=1;r<=Grate[2].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*Grate[2].Growth[r]/(coagulation_units*(float)Grate[2].NoEnsambles[r]));
		}
		fclose(taza);
		
		sprintf(distT,"%s/MetabolicR",contenedor);
		taza=fopen(distT, "w");
		for(r=1;r<=Grate[3].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*Grate[3].Growth[r]/(coagulation_units*(float)Grate[3].NoEnsambles[r]));
		}
		fclose(taza);
		
		//sprintf(distT,"%s/TotalResourceR",contenedor);
		//taza=fopen(distT, "w");
		//for(r=1;r<=Grate[4].i_max;r++){
			//fprintf(taza,"%d %f\n",r, Grate[4].Growth[r]/(coagulation_units*(float)Grate[4].NoEnsambles[r]));
		//}
		//fclose(taza);
	
FreeRate_log(&Grate[0]);
FreeRate_log(&Grate[1]);
FreeRate_log(&Grate[2]);
FreeRate_log(&Grate[3]);
FreeRate_log(&Grate[5]);

////

//GuardaRhoVsT_MP(contenedor,&MP_RhoVsT_1,NULL);
//LiberaMemoriaFloat2D_MP(&MP_RhoVsT_1);
LiberaMemoriaFloat1D_MP(&TamDist_1);

char thinning[50];
sprintf(thinning,"%s/thinning",contenedor);
file=fopen(thinning, "w");
fputs("# T meanSize meanDensity fisicalTime totalResourceR\n",file);

for(r=1;r<=T_max;r++){
fprintf(file,"%d %f %f %f %f\n",r, delta_s*(meanSize.array[r]/NoEnsambles), meanDensity.array[r]/NoEnsambles, time_map[r], Grate[4].Growth[r]/(coagulation_units*(float)Grate[4].NoEnsambles[r]) );
}
fclose(file);

FreeRate_log(&Grate[4]);
return;
}
