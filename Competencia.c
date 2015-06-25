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
#ifdef VIRTUAL_GRID			///< The algorithm should be define at top of file "libPP_6.1.h"
printf("VIRUAL_GRID=TRUE\n");		///< Prints the algorithm being used "Virtual grid of resources"
#endif
#ifdef SOI
printf("SOI=TRUE\n");		///< Prints the algorithm being used "Zone of Influence"
#endif

/** This section initialize the simulation parameteres */
runDescriptor run;		///< structure containing the simulation parameters
run.X=200; 			 	///< X lenght of the simulated field
run.Y=run.X;			///< Y lenght of the simulated field
run.grid_units=1.0; 	///< convertion factor between lenght units to virtual grid cell lenght.
run.size_units=1.0; 	///< convertion factor between lenght units and virtual tree size units.
run.T_max=3001;			///< Number of monte carlo sweeps in the simulation
run.NoEnsambles=200;	///< Number of simulations in the ensemble of simulations.
int const write_interval=100;	///< Number of monte carlo sweeps between writting state to disk.

//in this case the Area_units and Length_units are units of the imaginary resources grid.
run.Model.ResourcesScale = 1.0; ///< convertion factor between lenght units of resource grid and simulation length.
run.Model.competitionAsymetry = 80.0;	///< The asymetry of competition "p".
double const Area_units=(run.grid_units*run.Model.ResourcesScale)*(run.grid_units*run.Model.ResourcesScale);
double const Length_units=run.grid_units*run.Model.ResourcesScale;

double coagulation_units = 1.0;
#ifdef SOI
coagulation_units = 500.0;
run.Model.coagulation_factor = coagulation_units*1.0;
#else
run.Model.coagulation_factor = coagulation_units*1.0*3.1416;
#endif

run.Model.coagulation_radio_exp=1.0;  ///< Alometric relation: exponent relating steem radio and canopy radio 
run.Model.coagulation_exp=2*run.Model.coagulation_radio_exp; 
run.Model.Cr=1.0;		///< Alometric relation: 
run.Model.coagulation_radio_factor=1.0; ///< Factor between resource intake radio and canopy radio 
run.Model.metabolic_exp=2.6666; ///< Alometric relation: exponent relating dbh and biomass (8/3)
run.Model.Cm=0.4;
run.Model.metabolic_factor=(run.Model.Cm)/pow(2.0*run.Model.Cr,8.0/3.0);
run.Model.health_factor=0.0; ///< Fraction of allowed missed biomass before individual dies.
run.Model.Cg=0.5;
run.Model.growth_constant=(3.0*pow(run.Model.Cr,8.0/3.0))/(pow(2.0,1.0/3.0)*run.Model.Cg); ///<  needed resources per unit size increse.
run.Model.resource_rate=1.0/coagulation_units; 
run.initialMeanDistance=run.grid_units*1;  ///< mean distance between individuals (when random pattern is generated)
run.initialMinSeparation=1;	///< Minimun separation between individual in a random pattern (filter)
run.Model.birth_rate=0.0;	///< birth rate of new individuals
run.Model.RadioBirth=15;	///< radius of new offsprings from parent 
run.Model.dead_rate=0.0;	///< instrinsic dead rate of individuals
run.Model.intra_coagulation=0.0; ///< deprecated
Individual indv;
indv.species=1;			///< species type
indv.size_float=2.5;	///< biomass of new individuals floating point variable
indv.size=run.size_units*10.0*indv.size_float;
indv.radio_float=R(indv, (&modelo));  ///< dbh radio of individuals
indv.radio=indv.radio_float;	
indv.metabolism=0;	///< deprecated, used to improve time step performance
indv.health=0; 		///< deprecated, used to track individual healt
omp_set_num_threads(4);		///< sets the number of threads used by omp

/** internal initialization variables **/
float const separation=(float)run.initialMeanDistance;
int T_max = run.T_max;
int NoEnsambles=run.NoEnsambles;
run.Model.meta_needs=SetMetaNeeds(run.Model, Area_units); 
run.Model.R=SetR(run.Model, Length_units);
run.Model.M=SetM(run.Model, Area_units);
model modelo;
modelo = run.Model;
int NDX=(run.X*run.grid_units)+1;
int NDY=(run.Y*run.grid_units)+1;
double delta_s=1.0/(run.size_units);
///************************************************///

/** Initialization of sturctures that keep the resoults **/
	
Float1D_MP TamDist_1;		///< Instance of Float1D_MP that will store the size distribution at each steep
		InicializaFloat1D_MP(&TamDist_1, T_max+10);		///< initialization function of Float1D structures.
		
Float1D_MP CumulativeTamDist_1;		///< instance of Float1D_MP that can store the comulative size distribution
		InicializaFloat1D_MP(&CumulativeTamDist_1, T_max+10);
		
Float1D_MP MP_CorrelacionG;			///< instance of Float1D_MP that can store a pair correlation 
	InicializaFloat1D_MP(&MP_CorrelacionG, NDX);  ///<  initialization
	ResetFloat1D_MP(&MP_CorrelacionG);		///< we ensure it is full of zeros

char contenedor[150];					
	sprintf(contenedor,"DATOS_TAM");		///< Name of the directory where data is going to be stored
	CreaContenedor(contenedor,run);			///< Creates the directory of the data and write a file with some parameters of the simulation in it. 
	
Float1D_MP meanDensity;			///< Instance of Float1D_MP to keep track of the mean density evolution
	InicializaFloat1D_MP(&meanDensity, T_max+10);	///<  initialization
	
Float1D_MP meanSize;			///< Instance of Float1D_MP to keep track of the mean size evolution
	InicializaFloat1D_MP(&meanSize, T_max+10);		///< initialization
	
float time_map[T_max+10];	///< To map monte-carlo sweep with fisical time
time_map[0]=0.0;		///< initialization

Rate_log Grate[6];			///< Structure that can be passed to the monte carlo sweep to keep track of instantaneous rates. (used to mesoure whatever you need)
InitRate_log(&Grate[5],10);	///< Initialization
InitRate_log(&Grate[0],50);	///< Initialization

FILE *file;		///< handle the files to write data to disk

/******************************************************/

	#pragma omp parallel			///< begin of parallel code. This code will be done for each thread 
	{
		int num_threads = omp_get_num_threads();	///< number of threads running
		int id = omp_get_thread_num();				///< id of 'this' thread
		int MaxPar = NoEnsambles/num_threads;		///< Number of simulation that each thread will handle
		#pragma omp master							///< The remainer simulations will be handle by the master thread
		{
				MaxPar+= NoEnsambles - MaxPar * num_threads;	 ///< Adding the remainer simulations to the master thread
		}
		estado e[MaxPar];						///< Struct "estado" stores the state of each simulation (positions, sizes, time)
		int Par;				///< Dummy variable used for looping.
		for(Par=0;Par<MaxPar;Par++)				///< for each simulated field
		{
			AlojaMemoria(&e[Par], NDX, NDY);	///< allocate memory for the field (state of the simualtion: positions, sizes, time)
			ResetEstado(&e[Par]);		///< initialization 
			e[Par].units=run.grid_units;	///< each simulation could have it's own units (not really used)
			e[Par].size_units=run.size_units; ///< each simulation could have it's own units (not really used)
		}
		init_JKISS(); ///< Initialization of the random number generator "JKISS"
			
		int i,j;		///< dummy variables used for looping
		
 /** INITIAL POSITIONS of individuals */
			for(Par=0;Par<MaxPar;Par++)			///< For each fiel (simulation) generate the INITIAL STATE
			{
		/** Random initial position with specific mean and minimun distance **/
		//	GeneraEstadoAleatorioTamano(&e[Par], 1.0/(float)run.initialMeanDistance , indv);	///< Random initial specifing the mean 
		//	FilterMinDistance(&e[Par], run.initialMinSeparation);	///< filter of minimun distance
		/**********************************************************/
		/** Grid initial positions **/
				//for(i=separation;i<NDX;i+=30*separation)
				//{
					//for(j=separation;j<NDY;j+=30*separation)
					//{
						//if(i+separation < NDX && j<NDY){
							//InsertIndividualAt(&e[Par],i,j,indv,1);
							//InsertIndividualAt(&e[Par],i+separation,j,indv,1);
						//}
					//}
				//}
		/*********************************************************/		
		/** Random initial positions specifing number of individuals */
				while(e[Par].ON < 55)
				{
					i=I_JKISS(1,NDX);
					j=I_JKISS(1,NDY);
					InsertIndividualAt(&e[Par],i,j,indv,0);
				}
		/**************************************************************/
			
		setMaxMetabolic(&e[Par],&modelo);		/// internal use. sets the frecuency of the poisson events in the simulation. 
		}

		Float1D_MP TamDist;				///< Temporal array to store size distribution to take advantage of each core cache. 
		InicializaFloat1D_MP(&TamDist, T_max+100);
			
	//////////////////////////////Barrido Monte CARLO:
	/** Initialization of temporal internal variables of each thread to improve memory efficiency **/
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
		
		indv.species=1;
		int s1,s2,r2;
		int done=0;
		float Interaction,s1val,s2val,f, Area,Iss,instantDeadRate,cumDeadRate;
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
			#pragma omp barrier
				if((i-(i/write_interval)*write_interval)==1)    //Inicializa cada write_interval
				{
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
					
						if(i==3001){
							int z;
							float val,Dprev,Dnext,yval_f,yval_b;
							for(z=1;z<TamDist_1.i_max;z++)
							{
								yval_f=(TamDist_1.array[z])*(pow(z*TamDist_1.index_units,(5.0/8.0))*(8.0/(10.0*2.0*4.7*3.0)));
								yval_b=(TamDist_1.array[z-2])*(pow((z-2)*TamDist_1.index_units,(5.0/8.0))*(8.0/(10.0*2.0*4.7*3.0)));
								Dprev=yval_f - yval_b;
								yval_f=(TamDist_1.array[z+2])*(pow((z+2)*TamDist_1.index_units,(5.0/8.0))*(8.0/(10.0*2.0*4.7*3.0)));
								yval_b=(TamDist_1.array[z])*(pow(z*TamDist_1.index_units,(5.0/8.0))*(8.0/(10.0*2.0*4.7*3.0)));
								Dnext=yval_f - yval_b;
								if(Dprev > 0.0 && Dnext <=0.0 && yval_b >1.0)
								{
									val = 2.0*4.7*pow(z*TamDist_1.index_units,(3.0/8.0));
									printf("yval = %f , max in = %f\n", yval_b ,val);
								}
							}
						}
						
					//	sprintf(distT,"CumulativeDT_%d",i);
					//	GuardaFloat1D_MP(contenedor,distT,&CumulativeTamDist_1);
					}
					////cumulative and instant deadRate
					//#pragma omp single
					//{
						//InitRate_log(&Grate[1],rate[1].i_max);
					//}
					//SumRate_log(&rate[1], &Grate[1]);
					//#pragma omp barrier
					//#pragma omp master
					//{
						//sprintf(distT,"%s/cumDeadR_%d",contenedor,i);
						////sprintf(distT,"%s/instantDeadR_%d",contenedor,i);
						//file=fopen(distT, "w");
						////if(( Grate[1].i_max - Grate[0].i_max ) > 0)
					////	{
					////		ReallocRate_log(&Grate[0], Grate[1].i_max - Grate[0].i_max );
					////	}
						//for(j=1;j<=Grate[1].i_max;j++){
							//if(Grate[1].NoEnsambles[j]>0)
							//{
							////instantDeadRate = (Grate[1].Growth[j]/((float)Grate[1].NoEnsambles[j])) - (Grate[0].Growth[j]/((float)Grate[0].NoEnsambles[j]));
							//cumDeadRate = (Grate[1].Growth[j]/((float)Grate[1].NoEnsambles[j]));
							//fprintf(file,"%f %f\n",((float)j)*delta_s, cumDeadRate);
							//}
						////	Grate[0].Growth[j]=Grate[1].Growth[j];
						////	Grate[0].NoEnsambles[j]=Grate[1].NoEnsambles[j];
						//}
						//fclose(file);
						//FreeRate_log(&Grate[1]);
					//}		
					//competitive and mature deadRate
					//if(done == 0 && meanSize.array[i] > 100*NoEnsambles )
					//{
						//#pragma omp single
						//{
							//InitRate_log(&Grate[1],rate[1].i_max);
						//}
						//SumRate_log(&rate[1], &Grate[1]);
						//#pragma omp barrier
						//#pragma omp master
						//{
							//sprintf(distT,"%s/competitiveDeadR",contenedor);
							//file=fopen(distT, "w");		
							//if(( Grate[1].i_max - Grate[0].i_max ) > 0)
							//{
								//ReallocRate_log(&Grate[0], Grate[1].i_max - Grate[0].i_max );
							//}
							//for(j=1;j<=Grate[1].i_max;j++){
								//if(Grate[1].NoEnsambles[j]>0)
								//{
								//cumDeadRate = (Grate[1].Growth[j]/((float)Grate[1].NoEnsambles[j]));
								//fprintf(file,"%f %f\n",((float)j)*delta_s, cumDeadRate);
								//}
								//Grate[0].Growth[j]=Grate[1].Growth[j];
								//Grate[0].NoEnsambles[j]=Grate[1].NoEnsambles[j];
							//}
							//fclose(file);
							//FreeRate_log(&Grate[1]);	
						//}
						//done=1;	
					//}
					//if(i >= (T_max - write_interval) && done == 1 )
					//{
						//#pragma omp single
						//{
							//InitRate_log(&Grate[1],rate[1].i_max);
						//}
						//SumRate_log(&rate[1], &Grate[1]);
						//#pragma omp barrier
						//#pragma omp master
						//{
							//sprintf(distT,"%s/matureDeadR",contenedor);
							//file=fopen(distT, "w");		
							//for(j=1;j<=Grate[1].i_max;j++){
								//if(Grate[1].NoEnsambles[j]>0)
								//{
									//if(Grate[0].i_max >= j)
									//{
									//Grate[1].Growth[j]=Grate[1].Growth[j] - Grate[0].Growth[j];
									//Grate[1].NoEnsambles[j]=Grate[1].NoEnsambles[j]-Grate[0].NoEnsambles[j];
									//}
								//cumDeadRate = (Grate[1].Growth[j]/((float)Grate[1].NoEnsambles[j]));
								//fprintf(file,"%f %f\n",((float)j)*delta_s, cumDeadRate);
								//}
							//}
							//fclose(file);
							//FreeRate_log(&Grate[1]);	
						//}
						//done=2;	
					//}
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
	//LiberaMemoriaFloat2D_MP(&MP_Corr2D);
	//LiberaMemoriaFloat1D_MP(&MP_Correlacion);

	#pragma omp single
	{
		FreeRate_log(&Grate[0]);
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
			fprintf(taza,"%f %f\n",((float)r)*delta_s, Grate[1].Growth[r]/((double)Grate[1].NoEnsambles[r]));
		}
		fclose(taza);
		
		sprintf(distT,"%s/ResourceR",contenedor);
		taza=fopen(distT, "w");
		for(r=1;r<=Grate[2].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*Grate[2].Growth[r]/((float)Grate[2].NoEnsambles[r]));
		}
		fclose(taza);
		
		sprintf(distT,"%s/MetabolicR",contenedor);
		taza=fopen(distT, "w");
		for(r=1;r<=Grate[3].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*Grate[3].Growth[r]/((float)Grate[3].NoEnsambles[r]));
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
free(run.Model.meta_needs);
return;
}
