#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "FContact2D.h"
#include "EntSal2D.h"

main(){
int NDX=500;
float Lambda= 1.65;
float fracON = 0.01;
int NoPasos=2000;
int NDY=NDX;

estado e;
int Sitios = (NDY*NDX + 1);
sitio SO[Sitios];
int np;
float q;
char contenedor[100];
FILE *RhoVsT;

	AlojaMemoria(&e,NDX,NDY);
	GeneraEstadoAleatorio(&e, fracON);
	RellenaSO(&e,SO);

	sprintf(contenedor,"CP2.2_Lado@%d_Lambda@%1.4f_qIni@%1.2f",NDX,Lambda,fracON);
	CreaContenedor(contenedor);
	RhoVsT=AbreRhoVsTEn(contenedor);
//	GuardaEstadoEn(contenedor,&e);
	ActualizaRhoVsT(&e,RhoVsT);
	

	for (np=1;np<=NoPasos;np++)
	{
		BarridoMC3(&e,SO,Lambda);
	//	GuardaEstadoEn(contenedor,&e);
		ActualizaRhoVsT(&e,RhoVsT);

		if(e.ON==0 && np<(NoPasos-3)){np=NoPasos-2;}
	}
	fclose(RhoVsT);
}

