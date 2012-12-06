/*
Copyright 2012 Jorge Velazquez
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include "GNA.h"
#include "libPP_6.1.h"



especie *parametros=NULL;
static float Max_Metabolic=0.0;
 clock_t start, end, start2, end2;
static double cpu_time_used=0.0;
static double tiempo2=0.0;

Grupo GRUPO_INI = { 0, 0, 0, -1 }; //TIPO = 0 : Todos los tipos, s = 0 : Todos los tamanos, NEG = 0 : No negacion, Numero de elementos en la ultima vez que se proceso (-1 : no procesado todavia)


void AlojaMemoriaEspecie(int tipo)
{	
static int Max_especie=-1;
		if(Max_especie<tipo)
		{
			parametros = (especie *)realloc(parametros, (tipo + 1)*sizeof(especie));
			if ( parametros == NULL ) 
			{
				printf("Error en memoria de parametros!\n");
				exit(0);
			}else{
				Max_especie=tipo;
				//printf("Max_especie:%d\n",Max_especie);
			}
			//printf("Memoria alojada especie %d, de tamano:%d\n", tipo,(int)(tipo + 1)*sizeof(especie));
		}
return;
}

void SetBirth(float L,int tipo)
{
	AlojaMemoriaEspecie(tipo);
	parametros[tipo].Birth=L;
return;
}

void SetCoagulation(float e,int tipo)
{
	AlojaMemoriaEspecie(tipo);
	parametros[tipo].Coagulation=e;
return;
}

void SetCoagulationIntra(float e,int tipo)
{
	AlojaMemoriaEspecie(tipo);
	parametros[tipo].CoagulationIntra=e;
return;
}

void SetDead(float d,int tipo)
{
	AlojaMemoriaEspecie(tipo);
	parametros[tipo].Dead=d;
return;
}

void SetRadioBirth(int rb, int tipo)
{
	AlojaMemoriaEspecie(tipo);
	parametros[tipo].RadioBirth=rb;
return;
}

void SetRadioCoa(int rc, int tipo)
{
	AlojaMemoriaEspecie(tipo);
	parametros[tipo].RadioCoa=rc;
return;
}

void SetRadioCoaIntra(int rc, int tipo)
{
	AlojaMemoriaEspecie(tipo);
	parametros[tipo].RadioCoaIntra=rc;
return;
}

void EscalaTiempoMetabolico(int tipo)
{
float Metabolizmo=0.0;
	if(parametros==NULL)
	{
		printf("No hay parametros a normalizar\n");
		return;
	}else{
		Metabolizmo=parametros[tipo].Birth;
		Metabolizmo+=parametros[tipo].Coagulation;
		Metabolizmo+=parametros[tipo].CoagulationIntra;
		Metabolizmo+=parametros[tipo].Dead;
		if(Max_Metabolic<Metabolizmo)
		{
			Max_Metabolic=Metabolizmo;
		//	printf("escalado Metabolizmo a %f especie:%d\n",Max_Metabolic,tipo);
		}
	}
return;
}

void AlojaMemoria(estado *es, int NDX, int NDY)
{
int *imem;
int *INDmem;
sitio *SOmem;
int *TIPOmem;
int *AGEmem;
int col, ile;
int k;
col = NDY + 1;
ile = NDX + 1;

	imem = (int *)malloc(col * ile * sizeof(int));
	INDmem = (int *)malloc(col * ile * sizeof(int));
	SOmem = (sitio *)malloc(col * ile * sizeof(sitio));
	TIPOmem = (int *)malloc(col * ile * sizeof(int));
	AGEmem = (int *)malloc(col * ile * sizeof(int));
	if (imem == NULL || INDmem == NULL || SOmem == NULL || TIPOmem == NULL || AGEmem == NULL)
	{	
		puts("\nNo se pudo alojar memoria para el estado\n");
		if(imem == NULL){puts("Memoria para rejilla fallo\n");}
		if(INDmem == NULL){puts("Memoria para inices de lista fallo\n");}
		if(SOmem == NULL){puts("Memoria Lista de Sitios fallo\n");}
		if(TIPOmem == NULL){puts("Memoria TIPOS fallo\n");}
		exit(0);
	}
	memset(imem,0,col * ile * sizeof(int));
	memset(INDmem,0,col * ile * sizeof(int));
	memset(SOmem,0,col * ile * sizeof(sitio));
	memset(TIPOmem,0,col * ile * sizeof(int));
	memset(AGEmem,0,col * ile * sizeof(int));
	es->s = (int **)malloc(ile * sizeof(int *));
	es->INDICE = (int **)malloc(ile * sizeof(int *));
	es->SO = SOmem;
	es->TIPO = (int **)malloc(ile * sizeof(int *));
	es->AGE = (int **)malloc(ile * sizeof(int *));
	if ((*es).s == NULL || (*es).INDICE == NULL || (*es).TIPO == NULL || (*es).AGE == NULL)
	{	
		puts("\nNo se pudo alojar memoria para los apuntadores de la memoria\n");
		exit(0);
	}
	for (k = 0; k < ile; k++)
	{
		es->s[k] = imem + (k * col);
		es->INDICE[k] = INDmem + (k * col);
		es->TIPO[k] = TIPOmem + (k * col);
		es->AGE[k] = AGEmem + (k * col);
	}
	es->NDX = NDX;
	es->NDY = NDY;
	es->T = 0;
	es->ON = 0;
	es->Max_Metabolic=Max_Metabolic;
	es->Meta_T = 0.0;
return;
}

void LiberaMemoria(estado *es)
{
	free(es->s[0]);
	free(es->INDICE[0]);
	free(es->TIPO[0]);
	free(es->SO);
	free(es->s);
	free(es->INDICE);
	free(es->TIPO);
return;
}

void ResetEstado(estado *es)
{
int col = es->NDY + 1;
int ile = es->NDX + 1;
	memset(*(es->s),0,col * ile * sizeof(int));
	memset(*(es->INDICE),0,col * ile * sizeof(int));
	memset((es->SO),0,col * ile * sizeof(sitio));
	memset(*(es->TIPO),0,col * ile * sizeof(int));
	memset(*(es->AGE),0,col * ile * sizeof(int));
	es->ON = 0;
	es->T = 0;
	es->Max_Metabolic=Max_Metabolic;
	es->Meta_T = 0.0;
}

void GeneraEstadoAleatorio(estado *es, float frac, int tipo)
{
int NDX = es->NDX;
int NDY = es->NDY;
float fNDX = es->NDX;
float fNDY = es->NDY;
int col = es->NDY + 1;
int ile = es->NDX + 1;
int **s = es->s; 
int **Tip = es->TIPO;
sitio *SO = es->SO;
int **INDICE = es->INDICE;
long int xrand;
long int yrand;
int n=0;
int N, Libres;
int i,j;


N = (int)((fNDX * fNDY) * frac);

Libres=((NDX * NDY) - (es->ON))-N;
	
	if(frac==1.0)
	{
		(es->ON)=0;
		for(i=1;i<=NDX;i++)
		{
			for(j=1;j<=NDY;j++)
			{
			s[i][j]=1;
			Tip[i][j]=tipo;
			(es->ON)++;
			SO[(es->ON)].i = i;
			SO[(es->ON)].j = j;
			INDICE[i][j]=(es->ON);	
			}
		}
	}
	
	if(Libres>1)
	{
		while(n<N)
		{ 
			i = I_JKISS(1, NDX);
			j = I_JKISS(1, NDY);		
			if(s[i][j]<=0){
				s[i][j]=1;
				Tip[i][j]=tipo;
				n++;
				(es->ON)++;
				SO[(es->ON)].i = i;
				SO[(es->ON)].j = j;
				INDICE[i][j]=(es->ON);
			}
		}
	}else{
		for(i=1;i<=NDX;i++)
		{
			for(j=1;j<=NDY;j++)
			{
				if(s[i][j]<=0){
				s[i][j]=1;
				Tip[i][j]=tipo;
				(es->ON)++;
				SO[(es->ON)].i = i;
				SO[(es->ON)].j = j;
				INDICE[i][j]=(es->ON);
				}
			}
		}
	}	
	
return;
}

void ActualizaRyC(estado *es, int N, int campo)
{
float Rand; 
float pDead, pCreacion, pCoagulation1, C, Dead, Birth, Coagulation1, Nada,pCoa2;
sitio vecino;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i=es->SO[N].i;
int j=es->SO[N].j;
sitio *SO = es->SO;
int radioCre;
int radioCoa;
int ho = (es->TIPO[i][j] * 10) - 5;	//Nicho0 cambiar cuando cambie numero de especies distinto a 50 (solo en nicho es necesario)
	
	Rand = F_JKISS();
	
	Dead=parametros[es->TIPO[i][j]].Dead; // modelo neutral y J-C
	// rdead = ((float)(es->TIPO[i][j]-500)/5000.0)+0.05; //Modelo Lottery
	Birth=parametros[es->TIPO[i][j]].Birth;
//	Dead=Birth - (Birth - parametros[es->TIPO[i][j]].Dead) * exp(-pow(campo - ho,2.0)/20000.0);  // funcion nicho0 remaster
	Coagulation1=parametros[es->TIPO[i][j]].Coagulation;	
	pCoa2=parametros[es->TIPO[i][j]].CoagulationIntra/Max_Metabolic;
	pDead=Dead/Max_Metabolic;	//Asignar Max_Metabolic, si no hay division entre cero.
	pCreacion = Birth/Max_Metabolic;
	pCoagulation1 = Coagulation1/Max_Metabolic;

	if(Rand<=pDead) //aniquilacion
	{	
		s[i][j]=0;
		SO[N]=SO[(es->ON)];
		es->INDICE[SO[es->ON].i][SO[es->ON].j]=N;
		(es->ON)--;	
		es->TIPO[i][j]=0;
		
	}else{  //creation o coagulacion o nada
		if(Rand<=(pDead + pCreacion)) //creation
		{
			radioCre=parametros[es->TIPO[i][j]].RadioBirth;
			EligeUniforme(i,j,radioCre,&vecino);
			
			if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
			if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
			if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
			if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
			
			if(s[vecino.i][vecino.j]<=0)
			{
				s[vecino.i][vecino.j]=1;
				(es->ON)++;
				SO[(es->ON)]=vecino;
				es->INDICE[vecino.i][vecino.j]=(es->ON);
				es->TIPO[vecino.i][vecino.j]=es->TIPO[i][j];
			}
			
		}else{ //coagulacion1 o 2 o nada
			if(Rand<=(pDead + pCreacion + pCoagulation1))  //coagulacion
			{
				radioCoa=parametros[es->TIPO[i][j]].RadioCoa;
				 EligeUniforme(i,j,radioCoa,&vecino);
				if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
				if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
				if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
				if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
				
				if(s[vecino.i][vecino.j]>=1 && (es->TIPO[vecino.i][vecino.j])!=(es->TIPO[i][j]))     // mato a los que NO SON de mi propia especie 
				{
					s[vecino.i][vecino.j]=0;
					SO[(es->INDICE[vecino.i][vecino.j])]=SO[(es->ON)];
					es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[vecino.i][vecino.j]);
					(es->ON)--;
					es->TIPO[vecino.i][vecino.j]=0;
				}
			}else{  //coagulacion 2 o nada
				if(Rand<=(pDead + pCreacion + pCoagulation1 + pCoa2)) //coagulation Inter
				{
					radioCoa=parametros[es->TIPO[i][j]].RadioCoaIntra;
					 EligeUniforme(i,j,radioCoa,&vecino);
					if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
					if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
					if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
					if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
					
					if(s[vecino.i][vecino.j]>=1 && (es->TIPO[vecino.i][vecino.j])==(es->TIPO[i][j]))     //  Solo mato a mi propia especie 
					{
						s[vecino.i][vecino.j]=0;
						SO[(es->INDICE[vecino.i][vecino.j])]=SO[(es->ON)];
						es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[vecino.i][vecino.j]);
						(es->ON)--;
						es->TIPO[vecino.i][vecino.j]=0;
					}
				}	
			}
		}
	}
return;
}


void BarrMCcRyCamp(estado *es)
{
//	start2 = clock();  //clock comentar!!
int Indice,i,vtipo;
float DT=0.0;
int campo=0;
	
	while(DT<1.0){
		if((es->ON) != 0){
		
			 DT+=1.0/(es->ON); 
			Indice = I_JKISS(1,(es->ON));
			//campo = pow(250000-pow(es->SO[Indice].j - 500,2),0.5);  //nicho0 remaster campo 0:500
				campo = 0;  //neutral y J-C y Heromyopia
			ActualizaRyC(es,Indice,campo);	
		}else{
			DT=2.0;
		}
	}
	
	
/*	for(i=1;i<=1;i++)
	{
		vtipo=I_JKISS(0,999);
		 InsertaIndividuosAleatorio(es, 4, vtipo); //Valores para nicho
		// InsertaIndividuosAleatorio(es, 40, vtipo); //Valores lottery
		//InsertaIndividuosAleatorio(es, 4, vtipo); //Valores J-C
		//InsertaIndividuosAleatorio(es, 2, vtipo); //valores para neutral
	}
*/
(es->T)++;			

//     end2 = clock();			//clock comentar!!
  //   tiempo2 = ((double) (end2 - start2))/ CLOCKS_PER_SEC;;  //clock comentar!!	
  //   printf("Tiempo en elejir vecinos en una corrida:%f, tiempo corrida:%f\n",cpu_time_used,tiempo2);
   //  cpu_time_used=0.0;
		
return;
}


void EligeUniforme(int i,int j,int radio, sitio *vecino)
{
	start = clock();  //clock comentar!!

	int x,y,R2;
	int diam = 2 * radio;
	int r2 = radio * radio;
	
		do{
		x=I_JKISS(0,diam) - radio;
		y=I_JKISS(0,diam) - radio;

		R2= x * x + y * y;
		}while( R2 > r2 || R2==0);
		
		vecino->i = i + x;
		vecino->j = j + y;
		
	 end = clock();			//clock comentar!!
     cpu_time_used += ((double) (end - start))/ CLOCKS_PER_SEC;;  //clock comentar!!
	 
return; 	
}

void InsertaIndividuosAleatorio(estado *es, int N, int tipo)
{
int NDX = es->NDX;
int NDY = es->NDY;
float fNDX = es->NDX;
float fNDY = es->NDY;
int col = es->NDY + 1;
int ile = es->NDX + 1;
int **s = es->s; 
int **Tip = es->TIPO;
sitio *SO = es->SO;
int **INDICE = es->INDICE;
long int xrand;
long int yrand;
int n=0;
int Libres;
int i,j;


Libres=((NDX * NDY) - (es->ON))-N;
	
	
	if(Libres>1)
	{
		while(n<N)
		{ 
			i = I_JKISS(1, NDX);
			j = I_JKISS(1, NDY);		
			if(s[i][j]<=0){
				s[i][j]=1;
				Tip[i][j]=tipo;
				n++;
				(es->ON)++;
				SO[(es->ON)].i = i;
				SO[(es->ON)].j = j;
				INDICE[i][j]=(es->ON);
			}
		}
	}else{
		for(i=1;i<=NDX;i++)
		{
			for(j=1;j<=NDY;j++)
			{
				if(s[i][j]<=0){
				s[i][j]=1;
				Tip[i][j]=tipo;
				(es->ON)++;
				SO[(es->ON)].i = i;
				SO[(es->ON)].j = j;
				INDICE[i][j]=(es->ON);
				}
			}
		}
	}	
	

return;
}



int CuentaEspecie(estado *es, int tipo)
{
	int ON=es->ON;
	int n;
	int total=0;
	int **TIPO=es->TIPO;
	
	for(n=1;n<=ON;n++)
	{
		if(TIPO[es->SO[n].i][es->SO[n].j]==tipo)
		{
			total++;
		}
	}
return total;
}


void InsertaIndividuoEn(estado *es,int i,int j,int tipo)
{
int **s = es->s; 
int **TIPO = es->TIPO;
sitio *SO = es->SO;
int **INDICE = es->INDICE;
	
				if(s[i][j]<=0){
				s[i][j]=1;
				TIPO[i][j]=tipo;
				(es->ON)++;
				SO[(es->ON)].i = i;
				SO[(es->ON)].j = j;
				INDICE[i][j]=(es->ON);
				}
return;
}

void ActualizaRhoVsT_MP(estado *es,Float2D_MP *RhoVsT,Dist_MP *Dist,char Option)	
{
int T=es->T;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;
int **TIPO=es->TIPO;
sitio *SO=es->SO;
float rho,rho_specie;
int n, NoEspecies;

if(RhoVsT!=NULL)
{
	NoEspecies = (RhoVsT->j_max);
}else
{
	NoEspecies=0;
	if(Option=='s')
	{
		for(n=1;n<=ON;n++)
		{		
			if(NoEspecies<s[SO[n].i][SO[n].j])
			{
				NoEspecies=s[SO[n].i][SO[n].j];
			}		
		}
	}else{
		for(n=1;n<=ON;n++)
		{
			if(NoEspecies<TIPO[SO[n].i][SO[n].j])
			{
				NoEspecies=TIPO[SO[n].i][SO[n].j];
			}
		}
	}
}

int tot=NoEspecies + 1;

int rhoVec[tot];  
memset(rhoVec,0,tot * sizeof(int));
int NoPart;

rho=((float)ON)/((float)(NDX*NDY));

if(RhoVsT!=NULL)
{
#pragma omp atomic
RhoVsT->array[T][0]+=rho;	// Tipo 0 es el total
}

if(Dist!=NULL && Dist->T!=T)
{
	#pragma omp critical
	{
		memset(Dist->array,0, Dist->i_max * sizeof(int));
		Dist->T=T;
		Dist->NoEnsambles=0;
	}
}

	if(NoEspecies!=0)
	{
		if(Option=='s')  //opcion 's' significa size: densidad de individuos segun su tamano
		{
			for(n=1;n<=ON;n++)
			{
				if(s[SO[n].i][SO[n].j]<NoEspecies)
				{
					rhoVec[s[SO[n].i][SO[n].j]]++;
				}else{
					//printf("Hubo tamanos que no se registraron en RhoVsT: %d\n",s[SO[n].i][SO[n].j]);
				}
			}
		}else{
			for(n=1;n<=ON;n++)
			{
			rhoVec[TIPO[SO[n].i][SO[n].j]]++;			//Segmentation Fault si exite un tipo mas grande que el tamano de RhoVsT en j
			}
		}

		for(n=1;n<tot;n++)
		{
				if(rhoVec[n]!=0)
				{
					rho_specie=((float)rhoVec[n])/((float)(NDX*NDY));
					if(RhoVsT!=NULL)
					{
						#pragma omp atomic
						RhoVsT->array[T][n]+=rho_specie;
					}
					
					if(Dist!=NULL)
					{
							NoPart=(int)(rho_specie/Dist->TamParticion);
							#pragma omp atomic
							Dist->array[NoPart]++;
					}
				}
		}
		if(Dist!=NULL)
		{
			#pragma omp atomic
			Dist->NoEnsambles++;
		}
	}
	
return;
}
	
void IniciaMemoriaFloat2D_MP(Float2D_MP *ARRAY)
{
		int col = ARRAY->j_max + 1;
		int ile = ARRAY->i_max + 1;
		float *fmem;
		fmem = (float *)calloc(col * ile, sizeof(float));
		if(fmem==NULL)
		{
			puts("No se pudo alojar memoria para Float2D_MP\n");
			exit(0);
		}
		ARRAY->array=(float **)malloc(ile * sizeof(float *));
		if(ARRAY->array==NULL)
		{
			puts("No se pudo alojar memoria para punteros Float2D_MP\n");
			exit(0);
		}
		int k;
		for (k = 0; k < ile; k++)
		{
		ARRAY->array[k] = fmem + (k * col);
		}
		
return;		
}

void LiberaMemoriaFloat2D_MP(Float2D_MP *ARRAY)
{
	free(*(ARRAY->array));
	free(ARRAY->array);
return;
}

void IniciaMemoriaInt2D_MP(Int2D_MP *ARRAY)
{
		int col = ARRAY->j_max + 1;
		int ile = ARRAY->i_max + 1;
		int *imem;
		imem = (int *)calloc(col * ile, sizeof(int));
		if(imem==NULL)
		{
			puts("No se pudo alojar memoria para Int2D_MP\n");
			exit(0);
		}
		ARRAY->array=(int **)malloc(ile * sizeof(int *));
		if(ARRAY->array==NULL)
		{
			puts("No se pudo alojar memoria para punteros Int2D_MP\n");
			exit(0);
		}
		int k;
		for (k = 0; k < ile; k++)
		{
		ARRAY->array[k] = imem + (k * col);
		}
		
return;		
}

void IniciaMemoriaDist_MP(Dist_MP *Dist)
{
	float Dx = Dist->xFin - Dist->xIni;
	int tam = (int)(Dx/Dist->TamParticion) + 1;
		Dist->array = (int *)calloc(tam, sizeof(int));
		if(Dist->array==NULL)
		{
			puts("No se pudo alojar memoria para IntDist_MP\n");
			exit(0);
		}
		Dist->i_max=tam - 1;
		
return;
}

void ResetFloat2D_MP(Float2D_MP *ARRAY)
{
	int i,j;
	for(i=0;i<=ARRAY->i_max;i++)
	{
		for(j=0;j<=ARRAY->j_max;j++)
		{
			ARRAY->array[i][j]=0.0;
		}
	}
	ARRAY->NoEnsambles=0;
	ARRAY->T=0;
	
return;
}


void ResetFloat1D_MP(Float1D_MP *ARRAY)
{
	int i;
		for(i=0;i<=ARRAY->i_max;i++)
		{
			ARRAY->array[i]=0.0;
		}
		ARRAY->NoEnsambles=0;
		ARRAY->T=0;
		
return;
}

void ResetDist_MP(Dist_MP *Dist)
{
	int i;
	for(i=0;i<=Dist->i_max;i++)
	{
		Dist->array[i]=0;
	}
	Dist->T=0;
	Dist->NoEnsambles=0;

return;
}

void SumaFloat2D_MP(Float2D_MP *Origen, Float2D_MP *Destino)
{
int totj=(Destino->j_max);
int totx=Destino->i_max;
int T,tipo;
	for(T=0;T<=totx;T++)
	{
		for(tipo=0;tipo<=totj;tipo++)
		{	
		#pragma omp atomic
		Destino->array[T][tipo]+=Origen->array[T][tipo];		//Segmentation fault si el tamano de Origen y Destino no coinciden
		}
	}
	#pragma omp atomic
	Destino->NoEnsambles+=Origen->NoEnsambles;
	
	Destino->T=Origen->T;

return;	
}

void SumaDist_MP(Dist_MP *Origen, Dist_MP *Destino)
{
	
	int tot=Destino->i_max;
	int x;
		for(x=0;x<=tot;x++)
		{	
		#pragma omp atomic
		Destino->array[x]+=Origen->array[x];					//Segmentation fault si el tamano de Origen y Destino no coinciden
		}
	
	Destino->T=Origen->T;
	#pragma omp atomic
	Destino->NoEnsambles+=Origen->NoEnsambles;
	
return;
}

void InicializaFloat2D_MP(Float2D_MP *Objeto, int i_max, int j_max, int NoEnsambles)
{
		Objeto->i_max=i_max;
		Objeto->j_max=j_max;  //No especies
		Objeto->NoEnsambles=NoEnsambles;
		Objeto->T=0;
		IniciaMemoriaFloat2D_MP(Objeto);
return;
}

void InicializaDist_MP(Dist_MP *Objeto, float TamParticion, float xIni, float xFin)
{
	Objeto->TamParticion=TamParticion;
	Objeto->NoEnsambles=0;
	Objeto->xIni=xIni;
	Objeto->xFin=xFin;
	IniciaMemoriaDist_MP(Objeto);
return;
}

void GeneraEstadoAleatorioTamano(estado *es, float frac, int tipo, int tamano)
{
	int NDX = es->NDX;
int NDY = es->NDY;
float fNDX = es->NDX;
float fNDY = es->NDY;
int col = es->NDY + 1;
int ile = es->NDX + 1;
int **s = es->s; 
int **Tip = es->TIPO;
sitio *SO = es->SO;
int **INDICE = es->INDICE;
long int xrand;
long int yrand;
int n=0;
int N, Libres;
int i,j;




if(tamano<=0)
{
	for(i=1;i<=NDX;i++)
	{
		for(j=1;j<=NDY;j++)
		{
			if(s[i][j]<0)
			{
				s[i][j]=0;
			}
		}
	}
	
N = (int)((float)((NDX * NDY) - (es->ON)) * frac );	
Libres=((NDX * NDY) - (es->ON))-N;
}else{
	N = (int)((fNDX * fNDY) * frac);
	Libres=((NDX * NDY) - (es->ON))-N;	
}
	
	if(frac==1.0)
	{
		if(tamano>=1)
		{
			(es->ON)=0;
			for(i=1;i<=NDX;i++)
			{
				for(j=1;j<=NDY;j++)
				{
				s[i][j]=tamano;
				Tip[i][j]=tipo;
				(es->ON)++;
				SO[(es->ON)].i = i;
				SO[(es->ON)].j = j;
				INDICE[i][j]=(es->ON);	
				}
			}
		}else{
			for(i=1;i<=NDX;i++)
			{
				for(j=1;j<=NDY;j++)
				{
					if(s[i][j]<=0)
					{
					s[i][j]=tamano;
					Tip[i][j]=tipo;	
					}
				}
			}	
		}
	}
	
	if(Libres>1)
	{
		if(tamano>0)
		{
			while(n<N)
			{ 
				i = I_JKISS(1, NDX);
				j = I_JKISS(1, NDY);		
					if(s[i][j]<=0){
						s[i][j]=tamano;
						Tip[i][j]=tipo;
						n++;
						(es->ON)++;
						SO[(es->ON)].i = i;
						SO[(es->ON)].j = j;
						INDICE[i][j]=(es->ON);
					}
				
			}
		}else{
			while(n<N)
			{ 
				i = I_JKISS(1, NDX);
				j = I_JKISS(1, NDY);		
					if(s[i][j]<=0){
						s[i][j]=tamano;
						Tip[i][j]=tipo;
						n++;
					}
			}		
		}
	}else{
		if(tamano>0)
		{
			for(i=1;i<=NDX;i++)
			{
				for(j=1;j<=NDY;j++)
				{
					if(s[i][j]<=0){
					s[i][j]=tamano;
					Tip[i][j]=tipo;
					(es->ON)++;
					SO[(es->ON)].i = i;
					SO[(es->ON)].j = j;
					INDICE[i][j]=(es->ON);
					}
				}
			}
		}else{
			for(i=1;i<=NDX;i++)
			{
				for(j=1;j<=NDY;j++)
				{
					if(s[i][j]<=0){
					s[i][j]=tamano;
					Tip[i][j]=tipo;
					}
				}
			}
		}
	}	
	
	
return;
}

void ActualizaRyCTamano(estado *es, int N, int campo)
{
	float Rand; 
float pDead, pCreacion, pCoagulation1, C, Dead, Birth, Coagulation1, pCrecer;
sitio vecino;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i=es->SO[N].i;
int j=es->SO[N].j;
sitio *SO = es->SO;
int radioCre;
int radioCoa;
int ho = (es->TIPO[i][j] * 10) - 5;	//Nicho0 cambiar cuando cambie numero de especies distinto a 50 (solo en nicho es necesario)
	

//float bmin=14.0;
float NMax_Metabolic; 
	
	Rand = F_JKISS();
	
	
	Dead=parametros[es->TIPO[i][j]].Dead; // modelo neutral y J-C
	
	Birth=parametros[es->TIPO[i][j]].Birth;
	
//	pCoa2=parametros[es->TIPO[i][j]].CoagulationIntra/Max_Metabolic;
	pDead=Dead/es->Max_Metabolic;	//Asignar Max_Metabolic, si no hay division entre cero.
	pCreacion = Birth/es->Max_Metabolic;
	
	int max_tamano = 200;
	if(s[i][j]<=max_tamano)
	{
		Coagulation1 = (float)s[i][j];
	}else{
		Coagulation1 = (float)max_tamano;
	}
	
	//Coagulation1 = parametros[es->TIPO[i][j]].Coagulation*(max_tamano*(1.0-exp(-((float)s[i][j])/bmin)));
	//Coagulation1 = (parametros[es->TIPO[i][j]].Coagulation*((float)(s[i][j])));
	pCoagulation1 = Coagulation1/es->Max_Metabolic;
	
	if(Rand<=pDead) //aniquilacion
	{	
		s[i][j]=0;
		SO[N]=SO[(es->ON)];
		es->INDICE[SO[es->ON].i][SO[es->ON].j]=N;
		(es->ON)--;	
		es->TIPO[i][j]=0;
		
	}else{  //creation o coagulacion o nada
		if(Rand<=(pDead + pCreacion)) //creation
		{
			radioCre=parametros[es->TIPO[i][j]].RadioBirth;
			EligeUniforme(i,j,radioCre,&vecino);
			
			if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
			if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
			if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
			if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
			
			if(s[vecino.i][vecino.j]<=0)
			{
				s[vecino.i][vecino.j]=1;
				(es->ON)++;
				SO[(es->ON)]=vecino;
				es->INDICE[vecino.i][vecino.j]=(es->ON);
				es->TIPO[vecino.i][vecino.j]=es->TIPO[i][j];
			}
			
		}else{ //como o nada
			if(Rand<=(pDead + pCreacion + pCoagulation1))  //coagulacion
			{	
				radioCoa=parametros[es->TIPO[i][j]].RadioCoa;
				 EligeUniforme(i,j,radioCoa,&vecino);
				if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
				if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
				if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
				if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
				
				if(s[vecino.i][vecino.j]<=0)  //Si hay comida, como. 
				{
					s[vecino.i][vecino.j]=0;			
							s[i][j]++;			
							NMax_Metabolic = Birth + Dead + Coagulation1 + parametros[es->TIPO[i][j]].CoagulationIntra;
							if(es->Max_Metabolic < NMax_Metabolic)
							{
								es->Max_Metabolic = NMax_Metabolic;
							}
				}else{
					s[vecino.i][vecino.j]=0;
					SO[es->INDICE[vecino.i][vecino.j]]=SO[(es->ON)];
					es->INDICE[SO[es->ON].i][SO[es->ON].j]=es->INDICE[vecino.i][vecino.j];
					(es->ON)--;	
					es->TIPO[vecino.i][vecino.j]=0;
					es->AGE[vecino.i][vecino.j]=0;		
				}
			}
		}
	}
	
	
	
return;
}

void BarrMCcRyCampTamano(estado *es,float flujo_recursos,model *param)
{
int Indice,i,vtipo;
float DT=0.0;
float TMetabolicIni=es->Meta_T;
float LMax_Metabolic=es->Max_Metabolic;
float TEnteroAnterior;
float TMetabolicActual;
	
TEnteroAnterior = floor(es->Meta_T);
	
	while(DT<1.0){
		if((es->ON) != 0){
			TMetabolicActual= TMetabolicIni + DT/LMax_Metabolic;
			if(TMetabolicActual - TEnteroAnterior >=1.0)
			{
				GeneraEstadoAleatorioTamano(es, flujo_recursos, -1, -1);
				TEnteroAnterior += 1.0;
			}
			 DT+=1.0/(es->ON); 
			Indice = I_JKISS(1,(es->ON));
			//#pragma omp master
			//{
				//if(DT>=1.0){
			//printf("llega DT:%f\n",DT);
				//}
			//}
			ActualizaUniv(es, Indice, param);
			//ActualizaRyCTamano(es, Indice, 0);
		}else{
			DT=2.0;
		}
	}
	

(es->T)++;
es->Meta_T=TMetabolicActual;

return;
}

void ActualizaDistTamano_MP(estado *es, Float1D_MP *TamDist, char Opcion)
{
int T=es->T;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;
int **s = es->s;
sitio *SO=es->SO;
float prob_tam;
int n;
int max_tam=0;

		for(n=1;n<=ON;n++)
		{
			if(max_tam<s[SO[n].i][SO[n].j])
			{
				max_tam=s[SO[n].i][SO[n].j];
			}			
		}

int tot=max_tam + 1;

if(Opcion=='R')
{
	tot=(int)sqrt((double)max_tam) + 1; 
}

int rhoVec[tot];
memset(rhoVec,0,tot * sizeof(int));

#pragma omp critical
{
	if(TamDist->T!=T)
	{
		memset(TamDist->array,0, TamDist->i_max * sizeof(int));
		TamDist->T=T;
		TamDist->NoEnsambles=0;
	}
}
	if(Opcion=='R')
	{
		for(n=1;n<=ON;n++)
		{
			
			rhoVec[(int)sqrt((double)s[SO[n].i][SO[n].j])]++;			
		}
	}else{
		for(n=1;n<=ON;n++)
		{
			rhoVec[s[SO[n].i][SO[n].j]]++;			
		}
	}

		for(n=0;n<tot;n++)
		{
				if(rhoVec[n]!=0)
				{
					prob_tam=((float)rhoVec[n])/((float)(ON));
					#pragma omp atomic
					TamDist->array[n]+=prob_tam; //peligro semetation fault si 'n' es mas grande que el tamano de TamDist->array
				}
		}
		
		#pragma omp atomic
		TamDist->NoEnsambles++;
	
return;
}

void ActualizaRecursos_MP(estado *es,Float2D_MP *RhoVsT)
{
	int T=es->T;
	int n=RhoVsT->j_max;
	float recursos = 0.0;
	int i,j;
	int NDX=es->NDX;
	int NDY=es->NDY;
	
	for(i=1;i<=NDX;i++)
	{
		for(j=1;j<=NDY;j++)
		{
			if(es->s[i][j]<0)
			{
				recursos+=1.0;
			}
		}
	}
	
	recursos/=(float)(NDX*NDY);
	#pragma omp atomic
	RhoVsT->array[T][n]+=recursos;
	return;
}

void SetSpecie(int NoEspecie, float Birth, float Coagulation, float Dead, float RadioBirth, float RadioCoa)
{
	SetBirth(Birth,NoEspecie);
	SetCoagulation(Coagulation,NoEspecie);
	SetCoagulationIntra(0.0,NoEspecie);
	SetDead(Dead,NoEspecie);
	SetRadioBirth(RadioBirth,NoEspecie);
	SetRadioCoa(RadioCoa,NoEspecie);
	EscalaTiempoMetabolico(NoEspecie);
			
return;
}

void SetSpecie2(int NoEspecie, float Birth, float Coagulation,float CoagulationIntra, float Dead, float RadioBirth, float RadioCoa, float RadioCoaIntra )
{
		//printf("Especie %d inicializando\n",NoEspecie);
	SetBirth(Birth,NoEspecie);
	SetCoagulation(Coagulation,NoEspecie);
	SetCoagulationIntra(CoagulationIntra,NoEspecie);
	SetDead(Dead,NoEspecie);
	SetRadioBirth(RadioBirth,NoEspecie);
	SetRadioCoa(RadioCoa,NoEspecie);
	SetRadioCoaIntra(RadioCoaIntra,NoEspecie);
	EscalaTiempoMetabolico(NoEspecie);
		//	printf("Especie %d lista\n",NoEspecie);
return;
}


void SumaFloat1D_MP(Float1D_MP *Origen,Float1D_MP *Destino)
{
	int R_max=Origen->i_max;
	int i;
	
	for(i=1;i<=R_max;i++)
	{
		#pragma omp atomic
		Destino->array[i]+=Origen->array[i];
	}
	#pragma omp atomic
	Destino->NoEnsambles+=Origen->NoEnsambles;
	
	Destino->T=Origen->T;
	
return;
}

void InicializaFloat1D_MP(Float1D_MP *Objeto, int i_max)
{
		Objeto->NoEnsambles=0;
		Objeto->i_max=i_max;
		Objeto->T=0;

		int tam = i_max + 1;
		Objeto->array = (float *)calloc(tam, sizeof(float));
		if(Objeto->array==NULL)
		{
			puts("No se pudo alojar memoria para Float1D_MP\n");
			exit(1);
		}	
	
return;
}

void LiberaMemoriaFloat1D_MP(Float1D_MP *Objeto)
{
	free(Objeto->array);
return;
}

void CFFT(estado *es, Float2D_MP *correlacion)
{
	printf("Entra a fft\n");
fftw_complex *out;
fftw_plan p, plan2;
int NDX = es->NDX;
int NDY = es->NDY;
int **s = es->s;
double *in;
int i,j;
int nyh = ( NDY / 2 ) + 1;

#pragma omp critical
{
 in = fftw_alloc_real( sizeof ( double ) * NDX * NDY );
//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * NDY);


out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * nyh);


p = fftw_plan_dft_r2c_2d ( NDX, NDY, in, out, FFTW_ESTIMATE );
}
	printf("Termina primer plan\n");
 
 for ( i = 0; i < NDX; i++ )
  {
    for ( j = 0; j < NDY; j++ )
    {
     // in[i*NDY+j][0] = s[i+1][j+1];
      //in[i*NDY+j][1] = 0;
      in[i*NDY + j] = (double)s[i+1][j+1];
    }
   // printf("in 1comlumna: %f\n",in[i*NDY + 10]);
  }

	printf("Comienza FFT\n");
fftw_execute(p); /* repeat as needed */
	printf("Termina FFT\n");
 for ( i = 0; i < NDX; i++ )
  {
    for ( j = 0; j < nyh; j++ )
    {
      out[i*nyh + j][0] = out[i*nyh + j][0]*out[i*nyh + j][0] + out[i*nyh + j][1]*out[i*nyh + j][1];
      out[i*nyh + j][1] = 0;
    }
  //  printf("out trans 1 col: %f\n",out[i*nyh + 10][0]);
  }
		printf("2plan\n");
		#pragma omp critical
		{
plan2 = fftw_plan_dft_c2r_2d ( NDX, NDY, out, in, FFTW_ESTIMATE );
		}
	printf("2plan listo\n");
printf("Comienza FFT 2\n");
 fftw_execute ( plan2 );
 printf("Termina FFT 2\n");
 
 for ( i = 0; i < NDX; i++ )
  {
    for ( j = 0; j < NDY; j++ )
    {
		 in[i*NDY + j]/=(double)(NDX*NDY);
      correlacion->array[i][j]=in[i*NDY + j]/(double)(es->ON);
    }
  }
 
 #pragma omp critical
 {
fftw_destroy_plan(p);
fftw_destroy_plan(plan2);
}
fftw_free(in);
fftw_free(out);

return;
}

void CFFT_MP(estado *es, int NoEnsambles, Float2D_MP *correlacion)
{

fftw_complex *out;
fftw_plan p, plan2;
int NDX = es[0].NDX;
int NDY = es[0].NDY;
double *in;
int i,j,n;
int nyh = ( NDY / 2 ) + 1;

#pragma omp critical
{
in = fftw_alloc_real( sizeof ( double ) * NDX * NDY );
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * nyh);

p = fftw_plan_dft_r2c_2d ( NDX, NDY, in, out, FFTW_ESTIMATE );
plan2 = fftw_plan_dft_c2r_2d ( NDX, NDY, out, in, FFTW_ESTIMATE );
}

 for(n=0;n<NoEnsambles;n++)
 {
	 for ( i = 0; i < NDX; i++ )
	  {
		for ( j = 0; j < NDY; j++ )
		{
		  in[i*NDY + j] = (double)es[n].s[i+1][j+1];
		}
	   // printf("in 1comlumna: %f\n",in[i*NDY + 10]);
	  }

	fftw_execute(p); /* repeat as needed */
		
	 for ( i = 0; i < NDX; i++ )
	  {
		for ( j = 0; j < nyh; j++ )
		{
		  out[i*nyh + j][0] = out[i*nyh + j][0]*out[i*nyh + j][0] + out[i*nyh + j][1]*out[i*nyh + j][1];
		  out[i*nyh + j][1] = 0;
		}
	  }
			
	 fftw_execute ( plan2 );
	 
	 for ( i = 0; i < NDX; i++ )
	  {
		for ( j = 0; j < NDY; j++ )
		{
			 in[i*NDY + j]/=(double)(es[n].ON);
			 in[i*NDY + j]/=(double)(es[n].ON);
			 #pragma omp atomic
		  correlacion->array[i][j]+=in[i*NDY + j];
		}
	  }
	  #pragma omp atomic
	  correlacion->NoEnsambles++;
  }
 
 correlacion->T=es->T;
 
 #pragma omp critical
{
fftw_destroy_plan(p);
fftw_destroy_plan(plan2);
}
fftw_free(in);
fftw_free(out);

return;
}

void CompactaCorrelacion(Float2D_MP *corr2D, Float1D_MP *corrRadial)
{
	int radio;
	int r_max; 
	int MasDer;
	int RadioInterior2;
	int Radio2;
	int Contados,i,j;
	
	DoblaCorrelacion(corr2D);
	
	corrRadial->array[0]=corr2D->array[0][0];
	corrRadial->array[1]=(corr2D->array[1][0] + corr2D->array[0][1])/2.0;
	
	if((corr2D->j_max) < (corr2D->i_max))
	{
		r_max=corr2D->j_max;
	}else{
		r_max=corr2D->i_max;
	}
	
	for(radio=2;radio<=r_max/2;radio++)
	{
		Contados=0;
		RadioInterior2 = (radio-1) * (radio-1);
		Radio2 = radio * radio;
		MasDer=radio;
		i=MasDer;
		for(j=0;j<=i;j++)
		{ 
			do{
				if((i*i + j*j)<=Radio2)
				{
								//Octante 1
							corrRadial->array[radio]+=corr2D->array[i][j];
							Contados++;
					
					if(i!=j)
					{
								//Octante 2					
						corrRadial->array[radio]+=corr2D->array[j][i];
							Contados++;					
					}	
									
				}else{
					MasDer--;
				}
				i--;
			}while((i*i + j*j)>RadioInterior2);
			i=MasDer;
		}
		corrRadial->array[radio]/=(float)Contados;
	}
	corrRadial->NoEnsambles=corr2D->NoEnsambles;
	corrRadial->T=corr2D->T;
	corrRadial->i_max = r_max/2;
return;
}

void DoblaCorrelacion(Float2D_MP *corr2D)
{
	int i,j;
	
	int NDX=corr2D->i_max;
	int NDY=corr2D->j_max;
	int nyred = ( NDY / 2 ) + 1;
	int nxred = ( NDX / 2 ) + 1;
	
	corr2D->array[NDX][NDY]=corr2D->array[0][0];
	corr2D->array[0][NDY]=corr2D->array[0][0];
	corr2D->array[NDX][0]=corr2D->array[0][0];
	
	for(i=0;i<nxred;i++)
	{
		for(j=0;j<nyred;j++)
		{
			if(i==0)
			{
				corr2D->array[NDX][j]=corr2D->array[0][j];
				corr2D->array[NDX][NDY-j]=corr2D->array[0][NDY-j];
			}
			if(j==0)
			{
				corr2D->array[i][NDY]=corr2D->array[i][0];
				corr2D->array[NDX-i][NDY]=corr2D->array[NDX-i][0];
			}
				corr2D->array[i][j]=(corr2D->array[i][j] + corr2D->array[NDX-i][j] + corr2D->array[i][NDY-j] + corr2D->array[NDX-i][NDY-j])/4.0;
		}
	}
	
return;
}
	
	
void CFFT_Tipos_MP(estado *es, int NoEnsambles, Float2D_MP *correlacion,int TipoOrigen, int TipoDestino)
{
fftw_complex *out1;
fftw_complex *out2;
fftw_plan p1, p2, plan2;
int NDX = es[0].NDX;
int NDY = es[0].NDY;
double *in1;
double *in2;
double Re,Im;
int i,j,n,on1,on2;
int nyh = ( NDY / 2 ) + 1;

#pragma omp critical
{
in1 = fftw_alloc_real( sizeof ( double ) * NDX * NDY );
in2 = fftw_alloc_real( sizeof ( double ) * NDX * NDY );

out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * nyh);
out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * nyh);

p1 = fftw_plan_dft_r2c_2d ( NDX, NDY, in1, out1, FFTW_ESTIMATE );
p2 = fftw_plan_dft_r2c_2d ( NDX, NDY, in2, out2, FFTW_ESTIMATE );
plan2 = fftw_plan_dft_c2r_2d ( NDX, NDY, out1, in1, FFTW_ESTIMATE );
}

 for(n=0;n<NoEnsambles;n++)
 {
	 on1=0;
	 on2=0;
	 
	 for ( i = 0; i < NDX; i++ )
	  {
		for ( j = 0; j < NDY; j++ )
		{
			if(es[n].TIPO[i+1][j+1]==TipoOrigen)
			{
				in1[i*NDY + j] = (double)es[n].s[i+1][j+1];
				on1++;
			}else{
				in1[i*NDY + j] = 0.0;
			}
			if(TipoDestino < 0)
			{
				if(es[n].TIPO[i+1][j+1] != (-1*TipoDestino) && es[n].TIPO[i+1][j+1] != 0)
				{
					in2[i*NDY + j] = (double)es[n].s[i+1][j+1];
					on2++;
				}else{
					in2[i*NDY + j] = 0.0;
				}
			}else{	
				if(es[n].TIPO[i+1][j+1] == TipoDestino)
				{
					in2[i*NDY + j] = (double)es[n].s[i+1][j+1];
					on2++;
				}else{
					in2[i*NDY + j] = 0.0;
				}
			}
		}
	   // printf("in 1comlumna: %f\n",in[i*NDY + 10]);
	    
	  }

	if(on1>0 && on2>0)
	{
		
		fftw_execute(p1); 
		fftw_execute(p2);
		
		 for ( i = 0; i < NDX; i++ )
		  {
			for ( j = 0; j < nyh; j++ )
			{
				Re=out1[i*nyh + j][0]*out2[i*nyh + j][0] + out1[i*nyh + j][1]*out2[i*nyh + j][1];
				Im=out1[i*nyh + j][1]*out2[i*nyh + j][0] - out1[i*nyh + j][0]*out2[i*nyh + j][1];
			  out1[i*nyh + j][0] = Re;
			  out1[i*nyh + j][1] = Im;
			}
		  }
			
		 fftw_execute ( plan2 );
		 
		 for ( i = 0; i < NDX; i++ )
		  {
			for ( j = 0; j < NDY; j++ )
			{
				 in1[i*NDY + j]/=(double)on1;
				 in1[i*NDY + j]/=(double)on2;
				 #pragma omp atomic
					correlacion->array[i][j]+=in1[i*NDY + j];	
			}
		  }
		  
		  #pragma omp atomic
		  correlacion->NoEnsambles++;
		
		  correlacion->i_max=NDX;		 
		  correlacion->j_max=NDY;
		
	  }
	  
  }
 
 correlacion->T=es[0].T;
 
#pragma omp critical
{
fftw_destroy_plan(p1);
fftw_destroy_plan(p2);
fftw_destroy_plan(plan2);
fftw_free(in1);
fftw_free(in2);
fftw_free(out1);
fftw_free(out2);
}
return;
}


void CFFT_Univ_MP(estado *es, CorrDescriptor *Especifica, Float2D_MP *correlacion, Grupo *TipoOrigen, Grupo *TipoDestino)
{
fftw_complex *out1;
fftw_complex *out2;
fftw_plan p1, p2, plan2;
int NDX = es[0].NDX;
int NDY = es[0].NDY;
double *in1;
double *in2;
double Re,Im;
int i,j,n,on1,on2,on0, Muestra, ind;
int nyh = ( NDY / 2 ) + 1;
sitio SO[es[0].ON+1000];

int MeanSquare=Especifica->MeanSquare;
int NoEnsambles=Especifica->NoEnsambles;
int NoMuestras=Especifica->NoMuestras;
int MuestraIni;
int MuestraFin;
if (Especifica->Muestra > 0)	//Si se especifica Muestra, entonces solo se analiza esa muestra
{
	MuestraIni=Especifica->Muestra;
	MuestraFin=Especifica->Muestra;
}else{		// Valor de 0 (o negativo) en Especifica.Muestra analiza todas las muestras
	MuestraIni=1;
	MuestraFin=Especifica->NoMuestras;
}


#pragma omp critical
{
in1 = fftw_alloc_real( sizeof ( double ) * NDX * NDY );
in2 = fftw_alloc_real( sizeof ( double ) * NDX * NDY );

out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * nyh);
out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * nyh);

p1 = fftw_plan_dft_r2c_2d ( NDX, NDY, in1, out1, FFTW_ESTIMATE );
p2 = fftw_plan_dft_r2c_2d ( NDX, NDY, in2, out2, FFTW_ESTIMATE );
plan2 = fftw_plan_dft_c2r_2d ( NDX, NDY, out1, in1, FFTW_ESTIMATE );
}

 for(n=0;n<NoEnsambles;n++)
 {
	 on1=0;
	 on2=0;
	 on0=0;
	 
	 for ( i = 0; i < NDX; i++ )
	  {
		for ( j = 0; j < NDY; j++ )
		{
			if(TipoOrigen->TIPO!=0) //se sepecifico tipo
			{
				if(TipoOrigen->s!=0) //se especifico tamano
				{
					if(TipoOrigen->NEG!=1) //no es el negativo
					{
						if(es[n].TIPO[i+1][j+1]==TipoOrigen->TIPO && es[n].s[i+1][j+1]==TipoOrigen->s)
						{
							on0++;
							SO[on0].i=i;
							SO[on0].j=j;
						}
					}else{ //el negativo
						if(es[n].s[i+1][j+1]>0 && es[n].TIPO[i+1][j+1]!=TipoOrigen->TIPO && es[n].s[i+1][j+1]!=TipoOrigen->s)
						{
							on0++;
							SO[on0].i=i;
							SO[on0].j=j;
						}
					}
				}else{  //todos los tamanos
					if(TipoOrigen->NEG!=1) //no es el negativo
					{
						if(es[n].s[i+1][j+1]>0 && es[n].TIPO[i+1][j+1]==TipoOrigen->TIPO)
						{
							on0++;
							SO[on0].i=i;
							SO[on0].j=j;
						}
					}else{ //el negativo
						if(es[n].s[i+1][j+1]>0 && es[n].TIPO[i+1][j+1]!=TipoOrigen->TIPO)
						{
							on0++;
							SO[on0].i=i;
							SO[on0].j=j;
						}
					}
				}
			}else{ //todos los tipos
				if(TipoOrigen->s!=0) //se especifico tamano
				{
					if(TipoOrigen->NEG!=1) //no es el negativo
					{
						if(es[n].s[i+1][j+1]==TipoOrigen->s)
						{
							on0++;
							SO[on0].i=i;
							SO[on0].j=j;
						}
					}else{ //el negativo
						if(es[n].s[i+1][j+1]>0 && es[n].s[i+1][j+1]!=TipoOrigen->s)
						{
							on0++;
							SO[on0].i=i;
							SO[on0].j=j;
						}
					}
				}else{  //todos los tamanos
					if(TipoOrigen->NEG!=1) //no es el negativo
					{
						if(es[n].s[i+1][j+1]>0)
						{
							on0++;
							SO[on0].i=i;
							SO[on0].j=j;
						}
					}else{ //el negativo
						if(es[n].s[i+1][j+1]<1)
						{
							on0++;
							SO[on0].i=i;
							SO[on0].j=j;
						}
					}
				}
			}
			//--------------------
			if(TipoDestino->TIPO!=0) //se sepecifico tipo
			{
				if(TipoDestino->s!=0) //se especifico tamano
				{
					if(TipoDestino->NEG!=1) //no es el negativo
					{
						if(es[n].TIPO[i+1][j+1]==TipoDestino->TIPO && es[n].s[i+1][j+1]==TipoDestino->s)
						{
							in2[i*NDY + j] = 1.0;
							on2++;
						}else{
							in2[i*NDY + j] = 0.0;
						}
					}else{ //el negativo
						if(es[n].s[i+1][j+1]>0 && es[n].TIPO[i+1][j+1]!=TipoDestino->TIPO && es[n].s[i+1][j+1]!=TipoDestino->s)
						{
							in2[i*NDY + j] = 1.0;
							on2++;
						}else{
							in2[i*NDY + j] = 0.0;
						}
					}
				}else{  //todos los tamanos
					if(TipoDestino->NEG!=1) //no es el negativo
					{
						if(es[n].s[i+1][j+1]>0 && es[n].TIPO[i+1][j+1]==TipoDestino->TIPO)
						{
							in2[i*NDY + j] = 1.0;
							on2++;
						}else{
							in2[i*NDY + j] = 0.0;
						}
					}else{ //el negativo
						if(es[n].s[i+1][j+1]>0 && es[n].TIPO[i+1][j+1]!=TipoDestino->TIPO)
						{
							in2[i*NDY + j] = 1.0;
							on2++;
						}else{
							in2[i*NDY + j] = 0.0;
						}
					}
				}
			}else{ //todos los tipos
				if(TipoDestino->s!=0) //se especifico tamano
				{
					if(TipoDestino->NEG!=1) //no es el negativo
					{
						if(es[n].s[i+1][j+1]==TipoDestino->s)
						{
							in2[i*NDY + j] = 1.0;
							on2++;
						}else{
							in2[i*NDY + j] = 0.0;
						}
					}else{ //el negativo
						if(es[n].s[i+1][j+1]>0 && es[n].s[i+1][j+1]!=TipoDestino->s)
						{
							in2[i*NDY + j] = 1.0;
							on2++;
						}else{
							in2[i*NDY + j] = 0.0;
						}
					}
				}else{  //todos los tamanos
					if(TipoDestino->NEG!=1) //no es el negativo
					{
						if(es[n].s[i+1][j+1]>0)
						{
							in2[i*NDY + j] = 1.0;
							on2++;
						}else{
							in2[i*NDY + j] = 0.0;
						}
					}else{ //el negativo
						if(es[n].s[i+1][j+1]<1)
						{
							in2[i*NDY + j] = 1.0;
							on2++;
						}else{
							in2[i*NDY + j] = 0.0;
						}
					}
				}
			}
					
		}  
	  }
	  
	  
	  
	  if(MuestraIni != MuestraFin && MuestraFin>on0)
	  {
		  MuestraFin=on0;
		  printf("No maximo de muestras posibles = %d !\n",on0);
	  }
	  
	for(Muestra=MuestraIni;Muestra<=MuestraFin;Muestra++)
	{
		on1=0;
		for(i=0;i<(NDX*NDY);i++)
		{
			in1[i]=0.0;
		}
		for(ind=0;(ind*NoMuestras + Muestra)<=on0;ind++)
		{
			in1[SO[ind*NoMuestras + Muestra].i*NDY + SO[ind*NoMuestras + Muestra].j]=1.0;
			on1++;
		}

		if(on1>0 && on2>0)
		{
			
			fftw_execute(p1); 
			fftw_execute(p2);
			
			 for ( i = 0; i < NDX; i++ )
			  {
				for ( j = 0; j < nyh; j++ )
				{
					Re=out1[i*nyh + j][0]*out2[i*nyh + j][0] + out1[i*nyh + j][1]*out2[i*nyh + j][1];
					Im=out1[i*nyh + j][1]*out2[i*nyh + j][0] - out1[i*nyh + j][0]*out2[i*nyh + j][1];
				  out1[i*nyh + j][0] = Re;
				  out1[i*nyh + j][1] = Im;
				}
			  }
				
			 fftw_execute ( plan2 );
			 
			 for ( i = 0; i < NDX; i++ )
			  {
				for ( j = 0; j < NDY; j++ )
				{
					 in1[i*NDY + j]/=(double)on1;
					 in1[i*NDY + j]/=(double)on2;
					 #pragma omp atomic
						correlacion[0].array[i][j]+=in1[i*NDY + j];
						
					if(MeanSquare==1)
					{
						#pragma omp atomic
						correlacion[1].array[i][j]+=pow(in1[i*NDY + j],2.0);
					}
				}
			  }
		  
		  #pragma omp atomic
		  correlacion->NoEnsambles++;
		
		  correlacion->i_max=NDX;		 
		  correlacion->j_max=NDY;
		
		}
	}  
  }
 
 correlacion->T=es[0].T;
 TipoOrigen->on=on0;
 TipoDestino->on=on2;
 
 #pragma omp critical
{
fftw_destroy_plan(p1);
fftw_destroy_plan(p2);
fftw_destroy_plan(plan2);
fftw_free(in1);
fftw_free(in2);
fftw_free(out1);
fftw_free(out2);
}

return;
}

void ActualizaUniv(estado *es, int N, model *modelo)
{
float Rand; 
float pDead, pCreacion, pCoagulation1, Dead, Birth, CoagulationIntra, Coagulation1, pMetabolic;
sitio vecino;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i=es->SO[N].i;
int j=es->SO[N].j;
sitio *SO = es->SO;
int radioCre;
int radioCoa;
	
int max_tamano = 50;

float NMax_Metabolic; 
	
	Rand = F_JKISS();
	
	
	Dead=parametros[es->TIPO[i][j]].Dead; 
	Birth=parametros[es->TIPO[i][j]].Birth;
	CoagulationIntra=parametros[es->TIPO[i][j]].CoagulationIntra;

	if(s[i][j]<=20)
	{
	Coagulation1 = (float)s[i][j];
	}else{
		Coagulation1 = 0.5*(float)(s[i][j]-20) + 20.0;
	}
	
	pCreacion = Birth/es->Max_Metabolic;
	pDead=Dead/es->Max_Metabolic;
	pCoagulation1 = Coagulation1/es->Max_Metabolic;
	
	if(Rand<=pDead) //aniquilacion
	{	
		s[i][j]=0;
		SO[N]=SO[(es->ON)];
		es->INDICE[SO[es->ON].i][SO[es->ON].j]=N;
		(es->ON)--;	
		es->TIPO[i][j]=0;
		
	}else{  //creation o coagulacion o nada
		if(Rand<=(pDead + pCreacion)) //creation
		{
			radioCre=parametros[es->TIPO[i][j]].RadioBirth;
			EligeUniforme(i,j,radioCre,&vecino);
			
			if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
			if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
			if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
			if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
			
			if(s[vecino.i][vecino.j]<=0)
			{
				s[vecino.i][vecino.j]=1;
				(es->ON)++;
				SO[(es->ON)]=vecino;
				es->INDICE[vecino.i][vecino.j]=(es->ON);
				es->TIPO[vecino.i][vecino.j]=es->TIPO[i][j];
			}
			
		}else{ //como o nada
			if(Rand<=(pDead + pCreacion + pCoagulation1))  //coagulacion
			{	
				radioCoa=parametros[es->TIPO[i][j]].RadioCoa;
				 EligeUniforme(i,j,radioCoa,&vecino);
				if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
				if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
				if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
				if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
				
				//pMetabolic=(modelo->MetFact*pow((float)(s[i][j]),modelo->MetExp))/es->Max_Metabolic;
				//if(s[i][j]>20)
				//{
					//pMetabolic=1.5*(float)s[i][j]/es->Max_Metabolic;
				//}else{
					pMetabolic=0.5*(float)s[i][j]/es->Max_Metabolic;
				//}
				
				
				if(Rand<=(pDead + pCreacion + pMetabolic))  //Cuento las necesidades metabolicas.
				{
					es->AGE[i][j]++;
				}
				
				if(s[vecino.i][vecino.j]<0)  //Si hay comida, como. 
				{
					s[vecino.i][vecino.j]=0;
						if(es->AGE[i][j]<0)		//Si he llenado las necesidades de metabolizmo crezco
						{
							s[i][j]++;
							es->AGE[i][j]=0;
						
							NMax_Metabolic = Birth + Dead + Coagulation1 + CoagulationIntra;
							if(es->Max_Metabolic < NMax_Metabolic)
							{
								es->Max_Metabolic = NMax_Metabolic;
							}				
						}else{			//Si no he llenado necesidades de metabolizmo, como para ir llenando necesidades de metabolizmo.
							es->AGE[i][j]--;
						}
				}else{ //Si no hay comida, checo si corresponde morir.
					if(((modelo->ResurcesFact)*s[i][j]-es->AGE[i][j])<0) //Las reservas para satisfacer el metabolizmo, se proponen proporcionales al tamano. Si se acaban las reservas, muero.
					{
						s[i][j]=0;
						SO[N]=SO[(es->ON)];
						es->INDICE[SO[es->ON].i][SO[es->ON].j]=N;
						(es->ON)--;	
						es->TIPO[i][j]=0;
						es->AGE[i][j]=0;
					}
				}
			}
		}
	}	
	
return;
}

float LikelyHood(Float1D_MP *Origin, Float1D_MP *Experiment)
{
	float Total=0.0;
	float Sign;
	int i;
	int equilibrio=0;
	float FirstSign=0.0;
	float Coef=1.0;
	
	do{
		Sign=0;
		for(i=2;i<Experiment->i_max;i++)
		{
			Sign+=(Coef*Experiment->array[i])-(Origin->array[i]/Origin->NoEnsambles);
		}
		if(FirstSign==0.0)
		{
			FirstSign=Sign;
		}
		if(Sign*FirstSign<0.0 || Sign==0)
		{
			equilibrio=1;
		}else{
			if(Sign>0.0)
			{
				Coef-=0.01;
			}else{
				Coef+=0.01;
			}	
		}
	}while(equilibrio!=1);
	
		for(i=2;i<Experiment->i_max;i++)
		{
			Total+=fabs((Coef*Experiment->array[i])-(Origin->array[i]/Origin->NoEnsambles));
		}
	
	return Total;
}

void CargaExperiment(Float1D_MP *Experiment)
{
	Experiment->array[0]=0.0;
	Experiment->array[1]=0.0;
	Experiment->array[2]=0.0;
	Experiment->array[3]=1.0;
	Experiment->array[4]=2.0;
	Experiment->array[5]=3.0;
	Experiment->array[6]=6.0;
	Experiment->array[7]=8.0;
	Experiment->array[8]=10.0;
	Experiment->array[9]=13.0;
	Experiment->array[10]=15.0;
	Experiment->array[11]=16.0;
	Experiment->array[12]=17.0;
	Experiment->array[13]=19.0;
	Experiment->array[14]=21.0;
	Experiment->array[15]=20.0;
	Experiment->array[16]=17.0;
	Experiment->array[17]=17.0;
	Experiment->array[18]=18.0;
	Experiment->array[19]=18.0;
	Experiment->array[20]=17.0;
	Experiment->array[21]=19.0;
	Experiment->array[22]=21.0;
	Experiment->array[23]=23.0;
	Experiment->array[24]=21.0;
	Experiment->array[25]=19.0;
	Experiment->array[26]=18.0;
	Experiment->array[27]=17.0;
	Experiment->array[28]=16.0;
	Experiment->array[29]=15.0;
	Experiment->array[30]=13.0;
	Experiment->array[31]=9.0;
	Experiment->array[32]=6.0;
	Experiment->array[33]=5.0;
	Experiment->array[34]=4.0;
	Experiment->array[35]=3.0;
	Experiment->array[36]=2.0;
	Experiment->array[37]=1.5;
	Experiment->array[38]=1.0;
	Experiment->array[39]=0.7;
	Experiment->array[40]=0.3;
	Experiment->array[41]=0.0;
	
	int i;
	for(i=0;i<=41;i++)
	{
		Experiment->array[i]/=400.0;
	}
	
	Experiment->NoEnsambles=1;
return;
}


float Integra(Float1D_MP *Funcion, int inicial, int final)
{
	float Resultado=0.0;
	int i;
	for(i=inicial;i<=final;i++)
	{
		Resultado+=logf((Funcion->array[i])/((float)Funcion->NoEnsambles));
	}
	
return Resultado;
}

