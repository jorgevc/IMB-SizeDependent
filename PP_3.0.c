#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "libPP_2.0.h"
#include "EntSal2D.h"

main(){
int NDX=500;
float Lambda= 1.9;
float Epsilon=0.8;
float fracON = 0.09;
int NoPasos=500;
int NDY=NDX;

estado e;
int np;
char contenedor[100];
FILE *RhoVsT;
//	printf("direccionON antes de alojaMem=%d direccion e=%d\n",&e.ON,&e);

	AlojaMemoria(&e,NDX,NDY);
//	printf("e.ON = %d, e.NDX = %d direccionON=%d, direccion e=%d\n",e.ON,e.NDX,&e.ON,&e);
	GeneraEstadoAleatorio(&e, fracON);

	sprintf(contenedor,"PP_3.0_Lado@%d_Lambda@%1.4f_Epsilon@%1.4f_qIni@%1.2f",NDX,Lambda,Epsilon,fracON);
	CreaContenedor(contenedor);
	RhoVsT=AbreRhoVsTEn(contenedor);
	GuardaEstadoEn(contenedor,&e);
//	printf("e.ON = %d, e.NDX = %d direccionON=%d, direccion e=%d\n",e.ON,e.NDX,&e.ON,&e);
	ActualizaRhoVsT(&e,RhoVsT);

//	printf("e.ON = %d, e.NDX = %d direccionON=%d, direccion e=%d\n",e.ON,e.NDX,&e.ON,&e);
	

	for (np=1;np<=NoPasos;np++)
	{
		BarridoMC4(&e,Lambda,Epsilon);
		GuardaEstadoEn(contenedor,&e);
		ActualizaRhoVsT(&e,RhoVsT);

		if(e.ON==0 && np<(NoPasos-3)){np=NoPasos-2;}
	}
	fclose(RhoVsT);
}

