#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "FContact2D.h"
#include "EntSal2D.h"
#include "GNA.h"

main(){
int NDX=1000;
float Lambda= 1.65;
float DLambda = 0.01;
float fracON = 1;
int NoPasos=11000;
int NDY=NDX;

estado e;
int Sitios = (NDY*NDX + 1);
sitio SO[Sitios];

printf("Sitios:%d",Sitios);
int np;
int E;
int MAX_E = 20;
float rho_prom;
char archivo[100]="DATOS/";
unsigned int x,y,z,c;

char contenedor[100];
FILE *RhoVsL;

	sprintf(contenedor,"RvsL1.0_Lado@%d_T_Final@%d_qIni@%1.2f",NDX,NoPasos,fracON);
	CreaContenedor(contenedor);
	
	strcat(archivo,contenedor);
	strcat(archivo,"/RhoVsL");
	RhoVsL = fopen(archivo, "a");
		fputs("# Lambda   rho_st\n",RhoVsL);
		fclose(RhoVsL);

	AlojaMemoria(&e,NDX,NDY);

	while(Lambda<=1.65)
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
		RellenaSO(&e,SO);
	
			for (np=1;np<=NoPasos;np++)
			{
				BarridoMC3(&e,SO,Lambda);
			//	GuardaEstadoEn(contenedor,&e);
			//	ActualizaRhoVsT(&e,RhoVsT);
				if(e.ON==0 && np<(NoPasos-3)){np=NoPasos-2;}
			}
		rho_prom+=((float)e.ON)/((float)(NDX*NDY));
		printf("L= %f, Ensamble %d echo con semilla (%u , %u , %u , %u), rho_prom=%f \n",Lambda,E,x,y,z,c,rho_prom);
		}
		rho_prom = rho_prom/MAX_E;
	
		RhoVsL = fopen(archivo, "a");
		fprintf(RhoVsL,"%1.4f   %1.3f\n",Lambda,rho_prom);
		fclose(RhoVsL);
	Lambda+=DLambda;
	}
return;
}

