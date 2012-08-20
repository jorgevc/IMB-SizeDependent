#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "libPP_2.0.h"
#include "EntSal2D.h"
#include "GNA.h"

main(){
int NDX=1000;
float Lambda= 1.64;
float Epsilon = 0.0;
float DLambda = 0.01;
float DEpsilon = 0.01;
float fracON = 1.0;
int NoPasos=11000;
int NDY=NDX;

estado e;

int np; //indice numero de pasos
int E;  //indice numero de ensamble
int MAX_E = 7;
float rho_prom;
char archivo[100]="DATOS/";
char contenedor[100];
FILE *RhoVsL;


unsigned int x,y,z,c; //variables para la semilla

init_JKISS();  //inicializa con semilla aleatoria

	sprintf(contenedor,"RvsLxE_2.0_Lado@%d_T_Final@%d_qIni@%1.2f",NDX,NoPasos,fracON);
	CreaContenedor(contenedor);
	
	strcat(archivo,contenedor);
	strcat(archivo,"/RhoVsLxE.E3");
	RhoVsL = fopen(archivo, "a");
		fputs("# Lambda Epsilon rho_st\n",RhoVsL);
		fclose(RhoVsL);

	AlojaMemoria(&e,NDX,NDY);
	for(;Epsilon<=1.0;Epsilon+=DEpsilon)
	{
		for(;Lambda<=1.75;Lambda+=DLambda)
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
		
				for (np=1;np<=NoPasos;np++)
				{
					BarridoMC4(&e,Lambda,Epsilon);
				//	GuardaEstadoEn(contenedor,&e);
				//	ActualizaRhoVsT(&e,RhoVsT);
					if(e.ON==0 && np<(NoPasos-3)){np=NoPasos-2;}
				}
			rho_prom+=((float)e.ON)/((float)(NDX*NDY));
			printf("LxE= (%f,%f), E3.np= %d echo con semilla (%u , %u , %u , %u)\n",Lambda,Epsilon,E,x,y,z,c);
			}
			rho_prom = rho_prom/MAX_E;
			
			RhoVsL = fopen(archivo, "a");
			fprintf(RhoVsL,"%1.4f   %1.4f   %1.4f\n",Lambda, Epsilon, rho_prom);
			fclose(RhoVsL);
		}
	}
return;
}
