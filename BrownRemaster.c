#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "libPP_4.0.h"
#include "EntSalArb.h"
#include "GNA.h"

main(){
int NDX=1000;
float fracON = 1.0;
int NDY=NDX;
int T_max = 400;

float Birth1= 2.0;
float Coagulation1= 0.00002; 
float Dead1= 1.0;
int RadioBirth1= 10;
int RadioCoa1= 5;

SetBirth(3.0,0);
SetCoagulation(1.0,0);
SetDead(1.0,0);
//SetRadioBirth(RadioBirth1,NoEspecie);
//SetRadioCoa(RadioCoa1,NoEspecie);
EscalaTiempoMetabolico(0);
///////////////////////////////////// Estado INICIAL:
estado e;
//////////
//char contenedor2[150];
//sprintf(contenedor2,"BrRemasterTM_(B,D,RB)@(%f,%f,%d)_(NDX,Tmax,OnIni)@(%d,22000,5:5)",Birth,Dead,RadioBirth,NDX);
//CargaEstado(contenedor2, "T_3000", &e, NDX, NDY);
/////////

AlojaMemoria(&e, NDX, NDY);
ResetEstado(&e);

int NoEspecie;
for(NoEspecie=1;NoEspecie<=50;NoEspecie++)
{
SetBirth(Birth1,NoEspecie);
SetCoagulation(Coagulation1,NoEspecie);
SetDead(Dead1,NoEspecie);
SetRadioBirth(RadioBirth1,NoEspecie);
SetRadioCoa(RadioCoa1,NoEspecie);
// EscalaTiempoMetabolico(NoEspecie);
InsertaIndividuosAleatorio(&e,100,NoEspecie);
}


/////////////////////////////////////Termina Estado INICIAL
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS y guarda ESTADO INICIAL:
char contenedor[150];
sprintf(contenedor,"BrwRem_(B,D,C,RB,RC)@(%1.3f,%1.3f,%1.3f,%d,%d)_(NDX,Tmax)@(%d,%d)",Birth1,Dead1,Coagulation1,RadioBirth1,RadioCoa1,NDX,T_max);
CreaContenedor(contenedor);

FILE *datos;
datos=AbreRhoVsTEn(contenedor);
ActualizaRhoVsT(&e,datos,50);

//datos2=AbreNoSpeciesVsTEn(contenedor);

GuardaEstadoEn(contenedor, &e);

//ActualizaNoSpeciesVsT(datos2,GuardaTiposEn(contenedor, &e), e.T);

/////////////////////////////////////Termina Prepara CONTENEDOR para escribir DATOS y guarda ESTADO INICIAL
//////////////////////////////Barrido Monte CARLO:
int i;
for(i=1;i<=T_max;i++)
{
	BarrMCcRyCamp(&e);
	ActualizaRhoVsT(&e,datos,50);
	//GuardaEstadoEn(contenedor, &e);
	if((i-(i/20)*20)==0)    //Guarda Estado Cada 20 pasos
	{
		GuardaEstadoEn(contenedor, &e);
	}
}
GuardaEstadoEn(contenedor, &e);  //Guarda el ultimo estado de la corrida
fclose(datos);
///////////////////////////////Termina Monte CARLO

puts("Guardando Rango...\n");
GuardaTiposEn(contenedor, &e);
puts("Calculando Correlacion...\n");
GuardaCorrelacion(&e,1,250,contenedor);
puts("Calculando Correlacion entre 2 especies...\n");
GuardaCorrelacionTipo(&e,1, 250,1,25,contenedor);

return;
}
