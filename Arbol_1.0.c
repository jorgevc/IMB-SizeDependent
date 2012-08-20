#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "libPP_3.0.h"
#include "EntSalArb.h"
#include "GNA.h"

main(){
int NDX=500;
float fracON = 1;
int NDY=NDX;
int T_max = 5000;
float Birth= 1.63;
float Coagulation= 0.0; 
float Dead= 1.0;
int RadioBirth= 1;
int RadioCoa= 1;

char contenedor[100];

SetBirth(Birth);
SetCoagulation(Coagulation);
SetDead(Dead);
SetRadioBirth(RadioBirth);
SetRadioCoa(RadioCoa);

estado e;

AlojaMemoria(&e, NDX, NDY);
ResetEstado(&e);
GeneraEstadoAleatorio(&e,fracON, 1);

sprintf(contenedor,"A_Lado@%d_Tmax@%d_ONi@%1.1f_B-C-D-RB-RC@%1.3f-%1.3f-%1.3f-%d-%d",NDX,T_max,fracON,Birth,Coagulation,Dead,RadioBirth,RadioCoa);
CreaContenedor(contenedor);

GuardaEstadoEn(contenedor, &e);
int i;
for(i=1;i<=T_max;i++)
{
	BarrMCcRyCamp(&e);
	GuardaEstadoEn(contenedor, &e);
	printf("Tiempo: %d hecho\n",i);
}

return;
}
