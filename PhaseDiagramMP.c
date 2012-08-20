#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "libPP_5.0.h"
#include "EntSalArb_MP.h"
#include "GNA.h"
#include "ControlDinamico.h"

main(){	
//int clave=323;
//int *servicio = activa_escucha(clave);
//int *servicio;
//*servicio=0;
///////////////////////////Inicializa parametros de la simulacion
int NDX=200;
float fracON = 1.0;
int NDY=NDX;
int T_max = 600;
int NoEnsambles=100;

int CantidadEspecies=2;


float Birth1= 1.0;
float Coagulation1= 0.0; //Brown usa: 0.00002; 
float CoaIntra1= 0.2; //Modelo J-C 0.0008
float Dead1= 0.4;
int RadioBirth1= 10;
int RadioCoa1= 10;
int RadioCoaIntra1= 10;  //Modelo Heteromyopia 20


float Birth2= 1.0;
float Coagulation2= 0.0; //Brown usa: 0.00002;
float CoaIntra2=0.2;
float Dead2=0.3;
int RadioBirth2= 10;
int RadioCoa2= 10;
int RadioCoaIntra2= 10;


omp_set_num_threads(4);
////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:

Float2D_MP MP_RhoVsT_1;
int NoEspecies=CantidadEspecies;		
		InicializaFloat2D_MP(&MP_RhoVsT_1, T_max, CantidadEspecies, 0);
	

char contenedor[150];
	sprintf(contenedor,"PD(0.3:0.2)_(B,D,C,RB,RC,RCI)@(%1.3f,%1.3f:d,c:%1.4f,%d,%d,%d)_(NDX,Tmax)@(%d,%d)",Birth1,Dead1,Coagulation2,CoaIntra1,RadioBirth1,RadioCoa1,RadioCoaIntra1,NDX,T_max);
	//sprintf(contenedor,"BrwRemMPHM-2_(B,D,C,RB,RC)@(2.000,1.000,0.000,10,5)_(NDX,Tmax)@(1000,150)");
	CreaContenedor(contenedor);
	
char contenedorLec[150];	
	//sprintf(contenedorLec,"BrwRemNich0MP_(B,D,C,RB,RC)@(%1.3f,%1.3f,%1.3f,%d,%d)_(NDX,Tmax)@(%d,%d)",Birth1,Dead1,Coagulation1,RadioBirth1,RadioCoa1,NDX,T_max);
	sprintf(contenedorLec,contenedor);
/////////////////////////////////////Termina Prepara CONTENEDOR para escribir DATOS

char contenedorCompleto[300];

FILE *pD;
	char NombrePD[200];
	sprintf(NombrePD,"DATOS/%s/PD(d2,c1)",contenedor);
	pD=fopen(NombrePD, "a");
	fprintf(pD,"# D_2 C_1 Rho_1 Rho_2 Rho_tot\n"); 
	fclose(pD);

for(Dead2=0.40;Dead2<0.6;Dead2+=0.04)
{
	for(Coagulation1=0.0;Coagulation1<0.9;Coagulation1+=0.02)
	{
		
		SetSpecie2(1, Birth1, Coagulation1, CoaIntra1, Dead1, RadioBirth1, RadioCoa1, RadioCoaIntra1);
		SetSpecie2(2, Birth2, Coagulation2, CoaIntra2, Dead2, RadioBirth2, RadioCoa2, RadioCoaIntra2);
		
		sprintf(contenedorCompleto,"%s/(D2,C1)@(%1.3f,%1.3f)",contenedor,Dead2,Coagulation1);
		CreaContenedor(contenedorCompleto);
		
		ResetFloat2D_MP(&MP_RhoVsT_1);
		///////////////////////////////////// INICIA PARALLEL

		#pragma omp parallel			///////Estado INICIAL:
		{
			init_JKISS(); //Inicializa la semilla de cada proceso.
			
			int num_threads = omp_get_num_threads();
			int id = omp_get_thread_num();
			int MaxPar = NoEnsambles/num_threads;
			#pragma omp master
			{
					MaxPar+= NoEnsambles - MaxPar * num_threads;	 
			}
			estado e[MaxPar];
			
			//MaxPar=CargaEstado_MP(contenedorLec,"T_001",e,NDX,NDY,id,MaxPar);
			
			int Par;
			for(Par=0;Par<MaxPar;Par++)
			{
				AlojaMemoria(&e[Par], NDX, NDY);
				ResetEstado(&e[Par]); 	
			}

			
				for(Par=0;Par<MaxPar;Par++)
				{
				//InsertaIndividuosAleatorio(&e[Par],100,NoEspecie);
				GeneraEstadoAleatorio(&e[Par], 0.3, 1);
				GeneraEstadoAleatorio(&e[Par], 0.2, 2);
				}


			
		/////////////////////////////////////Termina Estado INICIAL
		//////////////////////////Prepara Contenedor en Memoria de cada proceso Para Mejorar rendimiento (optimizar el uso de cache de cada procesador)

			Float2D_MP MP_RhoVsT;	
					InicializaFloat2D_MP(&MP_RhoVsT, T_max, NoEspecies, MaxPar);
							
		/////////////////////////////////Termina prepara Contenedor MEMORIA de cada PROCESO

				for(Par=0;Par<MaxPar;Par++)
				{
					ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);	
				}
				
		//////////////////////////////Barrido Monte CARLO:
			int i;
			for(i=0;i<T_max;i++)
			{
				for(Par=0;Par<MaxPar;Par++)
				{
					//if((i-(i/10)*10)==0) // Guarda Estado cada 10 pasos
					//{
					//GuardaEstadoEn_MP(contenedor,&e[Par],id,Par);
					//}
					BarrMCcRyCamp(&e[Par]);
					ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
				}
				
				//if(*servicio==1)    //Servicio de control Dinamico
				//{
					//printf("entro a servicio\n");
					//Servicio(i,contenedor,&MP_RhoVsT, &MP_RhoVsT_1, NULL, NULL, e, MaxPar);
					//SalidaCD(&i,T_max);
				//} 
			
					if((i-(i/500)*500)==499)    //Inicializa cada 500 pasos
					{
						init_JKISS();
						printf("Se ha reinicializado JKISS (c/500pasos) paso:%d, threath:%d \n",i,id);
					}

			}
			
		//////////////////////////////Termina Monte CARLO
			
			SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
			PD_GuardaEstadoEn_MP(contenedorCompleto, e, id, MaxPar);
		

		}	////////////////////////////////////////////////////////////////////TERMINA PARALLEL


		GuardaRhoVsT_MP(contenedorCompleto,&MP_RhoVsT_1,NULL);
		
		

			pD=fopen(NombrePD, "a");
			fprintf(pD,"%f %f %f %f %f\n",Dead2,Coagulation1,MP_RhoVsT_1.array[T_max][1] , MP_RhoVsT_1.array[T_max][2], MP_RhoVsT_1.array[T_max][0]); 
			fclose(pD);
		
		
	}

}


//cierra_escucha();

return;
}
