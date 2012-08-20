#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "libPP_3.0.h"
#include "EntSalArb.h"
#include "GNA.h"

main(){
int NDX=1000;
float fracON = 0.0001;
int NDY=NDX;
int T_max = 1000;
float Birth= 0.4;
float Coagulation= 0.00004; 
float Dead= 0.2;
int RadioBirth= 10;
int RadioCoa= 5;

char contenedor[100];
FILE *datos2;

SetBirth(Birth);
SetCoagulation(Coagulation);
SetDead(Dead);
SetRadioBirth(RadioBirth);
SetRadioCoa(RadioCoa);

estado e;

sprintf(contenedor,"BrownNichoCorr");
CargaEstado(contenedor, "T_999", &e);
/*
AlojaMemoria(&e, NDX, NDY);
ResetEstado(&e);

int tipo,vtipo;
for(tipo=0;tipo<50;tipo++)
{
vtipo=I_JKISS(0,999);
InsertaIndividuosAleatorio(&e,50, vtipo);
}

sprintf(contenedor,"BrownNichoCorr");
//CreaContenedor(contenedor);

//datos=AbreRhoVsTEn(contenedor);
//datos2=AbreNoSpeciesVsTEn(contenedor);
//ActualizaRhoVsT(&e,datos);
GuardaEstadoEn(contenedor, &e);

//ActualizaNoSpeciesVsT(datos2,GuardaTiposEn(contenedor, &e), e.T);

int i;
for(i=1;i<=T_max;i++)
{
	BarrMCcRyCamp(&e);
//	ActualizaRhoVsT(&e,datos);
	GuardaEstadoEn(contenedor, &e);
	//GuardaTiposEn(contenedor, &e);
	//ActualizaNoSpeciesVsT(datos2,GuardaTiposEn(contenedor, &e), e.T);
	//printf("Tiempo: %d hecho\n",i);
	
}

//GuardaEstadoEn(contenedor, &e);
puts("calculando correlacion ...\n");
//GuardaCorrelacion(&e,1,250,contenedor);
*/

 GuardaCorrelacionTipo(&e,1, 250,804,713,contenedor);
//fclose(datos);


return;
}
