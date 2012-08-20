#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "libPP_2.0.h"
#include "EntSal2D.h"
#include "GNA.h"

main(){
int NDX=1000;
float LamIni= 1.648;
float EpsIni = 0.0;
float DLambda = 0.001;
float DEpsilon = 0.001;
float fracON = 1.0;
int NDY=NDX;


int MAX_E = 4;
int NoPasos=11000;
char archivo[100]="DATOS/";
char contenedor[100];
float rho_st=0.0;
float RvT[(NoPasos+2)];
int i; //inicializa RvT
 
 for(i=0;i<=(NoPasos+1);i++)
 {
 	RvT[i]=0.0;
 }
 

FILE *RhoVsL;

sprintf(contenedor,"RvsLxE_2.1_Lado@%d_T_Final@%d_qIni@%1.2f",NDX,NoPasos,fracON);
	CreaContenedor(contenedor);
	
	strcat(archivo,contenedor);
	strcat(archivo,"/RhoVsLxE.MP");
	RhoVsL = fopen(archivo, "a");
		fputs("# Lambda Epsilon rho_st\n",RhoVsL);
		fclose(RhoVsL);

#pragma omp parallel
  {
	estado e;
	int np; //indice numero de pasos
	int E;  //indice numero de ensamble
	float rho_prom, rho_t, Lambda,Epsilon;
	unsigned int x,y,z,c; //variables para la semilla
	int th_id = omp_get_thread_num();
	int no_th = omp_get_num_threads(); 
	FILE *RVL;
	FILE *RhoVsT;
	int T;
	char ArchRvT[100]="DATOS/";
	char nombreRvT[60];
	char nombreTemp[100];
	strcat(ArchRvT,contenedor);
	strcat(ArchRvT,"/");
	
		init_JKISS();  //inicializa con semilla aleatoria

		AlojaMemoria(&e,NDX,NDY);
		
	for(Epsilon=EpsIni;Epsilon<=0.1;Epsilon+=DEpsilon)
	{
		for(Lambda=LamIni;Lambda<=1.69;Lambda+=DLambda)
		{
			rho_prom=0.0;

			for(E=1;E<=MAX_E;E++)
			{
			x=JKISS();
			y=JKISS();
			z=JKISS();
			c=JKISS();
			Seed_JKISS(x,y,z,c);
			GeneraEstadoAleatorio(&e, fracON);
		
				for (np=0;np<NoPasos;np++)
				{
					rho_t=((float)e.ON)/((float)(e.NDX*e.NDY));
	                    rho_t=rho_t/(no_th*MAX_E);
	              
	                    #pragma omp atomic
					RvT[e.T]+=rho_t;
					
					//
					BarridoMC4(&e,Lambda,Epsilon);
				//	GuardaEstadoEn(contenedor,&e);
					if(e.ON==0 && np<(NoPasos-3)){np=NoPasos-2;}
				}			
					rho_t=((float)e.ON)/((float)(e.NDX*e.NDY));
	                    rho_t=rho_t/(no_th*MAX_E);
	                    #pragma omp atomic
					RvT[e.T]+=rho_t;	
					
			rho_prom+=((float)e.ON)/((float)(e.NDX*e.NDY));
			//printf("LxE= (%f,%f), %d.Ensamble= %d semilla (%u , %u , %u , %u)\n",Lambda,Epsilon,th_id,E,x,y,z,c);
			//printf("rho_prom = %d E=%d ID=%d\n",e.ON,E,th_id);
			}
			#pragma omp atomic 
			rho_st+=rho_prom;
			
			#pragma omp barrier
			
			//printf("Escribiendo datos... LxE= (%f,%f), %d.Ensamble= %d\n",Lambda,Epsilon,th_id,E);
			 
			 #pragma omp master
			 {
			 	rho_st = rho_st/(no_th*MAX_E);
		 		RVL = fopen(archivo, "a");
				fprintf(RVL,"%1.4f   %1.4f   %1.4f\n",Lambda, Epsilon, rho_st);
				fclose(RVL);
				rho_st = 0.0;
				//
				strcpy(nombreTemp,ArchRvT);
				sprintf(nombreRvT,"RvT.MP_Lambda@%1.4f_Epsilon@%1.4f",Lambda, Epsilon);
				strcat(nombreTemp,nombreRvT);
				RhoVsT = fopen(nombreTemp,"w");
				fputs("# T   rho \n",RhoVsT);
				for(T=0;T<=NoPasos;T++)
				{
					fprintf(RhoVsT,"%d   %1.4f\n",T,RvT[T]);
					RvT[T]=0.0;
				}
				fclose(RhoVsT);
			 }
			 
			 #pragma omp barrier
			
		}
	}
  }
return;
}
