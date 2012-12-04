/*
Copyright 2012 Jorge Velazquez
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "GNA.h"
#include "libPP_4.0.h"


especie *parametros=NULL;
static float Max_Metabolic=0.0;
 clock_t start, end, start2, end2;
static double cpu_time_used=0.0;
static double tiempo2=0.0;

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
			}
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

void EscalaTiempoMetabolico(int tipo)
{
float Metabolizmo;
	if(parametros==NULL)
	{
		puts("No hay parametros a normalizar\n");
		return;
	}else{
		Metabolizmo=parametros[tipo].Birth;
		Metabolizmo+=parametros[tipo].Coagulation;
		Metabolizmo+=parametros[tipo].Dead;
		if(Max_Metabolic<Metabolizmo)
		{
			Max_Metabolic=Metabolizmo;
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
int col, ile;
int k;
col = NDY + 1;
ile = NDX + 1;

	imem = (int *)malloc(col * ile * sizeof(int));
	INDmem = (int *)malloc(col * ile * sizeof(int));
	SOmem = (sitio *)malloc(col * ile * sizeof(sitio));
	TIPOmem = (int *)malloc(col * ile * sizeof(int));
	if (imem == NULL || INDmem == NULL || SOmem == NULL || TIPOmem == NULL)
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
	es->s = (int **)malloc(ile * sizeof(int *));
	es->INDICE = (int **)malloc(ile * sizeof(int *));
	es->SO = SOmem;
	es->TIPO = (int **)malloc(ile * sizeof(int *));
	if ((*es).s == NULL || (*es).INDICE == NULL || (*es).TIPO == NULL)
	{	
		puts("\nNo se pudo alojar memoria para los apuntadores de la memoria\n");
		exit(0);
	}
	for (k = 0; k < ile; k++)
	{
		es->s[k] = imem + (k * col);
		es->INDICE[k] = INDmem + (k * col);
		es->TIPO[k] = TIPOmem + (k * col);
	}
	es->NDX = NDX;
	es->NDY = NDY;
	es->T = 0;
	es->ON = 0;
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
	es->ON = 0;
	es->T = 0;
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
			if(s[i][j]==0){
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
				if(s[i][j]==0){
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
float pDead, pCreacion, pCoagulation, C, Dead, Birth, Coagulation, Nada;
sitio vecino;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i=es->SO[N].i;
int j=es->SO[N].j;
sitio *SO = es->SO;
int radioCre;
int radioCoa;
int ho = (es->TIPO[i][j] * 10) - 5;	//Nicho0 cambiar cuando cambie numero de especies distinto a 50
	
	Rand = F_JKISS();
	
	Dead=parametros[es->TIPO[i][j]].Dead; // modelo neutral y J-C
	// rdead = ((float)(es->TIPO[i][j]-500)/5000.0)+0.05; //Modelo Lottery
	Birth=parametros[es->TIPO[i][j]].Birth;
//	Dead=Birth - (Birth - parametros[es->TIPO[i][j]].Dead) * exp(-pow(campo - ho,2.0)/20000.0);  // funcion nicho0 remaster
	Coagulation=parametros[es->TIPO[i][j]].Coagulation;	
	pDead=Dead/Max_Metabolic;	//Asignar Max_Metabolic, si no hay division entre cero.
	pCreacion = Birth/Max_Metabolic;
	pCoagulation = Coagulation/Max_Metabolic;

	if(Rand<=pDead) //aniquilacion
	{	
		s[i][j]=0;
		SO[N]=SO[(es->ON)];
		es->INDICE[SO[es->ON].i][SO[es->ON].j]=N;
		(es->ON)--;	
		es->TIPO[i][j]=0;
		
	}else{  //creation o coagulacion
		if(Rand<=(pDead + pCreacion)) //creation
		{
			radioCre=parametros[es->TIPO[i][j]].RadioBirth;
			EligeUniforme(i,j,radioCre,&vecino);
			
			if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
			if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
			if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
			if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
			
			if(s[vecino.i][vecino.j]==0)
			{
				s[vecino.i][vecino.j]=1;
				(es->ON)++;
				SO[(es->ON)]=vecino;
				es->INDICE[vecino.i][vecino.j]=(es->ON);
				es->TIPO[vecino.i][vecino.j]=es->TIPO[i][j];
			}
			
		}else{ 
			if(Rand<=(pDead + pCreacion + pCoagulation))  //coagulacion
			{
				radioCoa=parametros[es->TIPO[i][j]].RadioCoa;
				 EligeUniforme(i,j,radioCoa,&vecino);
				if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
				if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
				if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
				if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
				
				if(s[vecino.i][vecino.j]==1)//  && (es->TIPO[vecino.i][vecino.j])!=(es->TIPO[i][j]))       No como a los de mi especie
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
return;
}


void BarrMCcRyCamp(estado *es)
{
//	start2 = clock();  //clock comentar!!
int Indice,i,vtipo;
float DT=0.0;
int campo;
	
	while(DT<1.0){
		if((es->ON) != 0){
		
			 DT+=1.0/(es->ON); 
			Indice = I_JKISS(1,(es->ON));
			//campo = pow(250000-pow(es->SO[Indice].j - 500,2),0.5);  //nicho0 remaster campo 0:500
				campo = 0;  //neutral y J-C
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
			if(s[i][j]==0){
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
				if(s[i][j]==0){
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

float OnPromRadio(estado *es, int radio)
{
	int x,y,R2,w,n;
	int diam = 2 * radio;
	int r2 = radio * radio;
	int A = diam * diam;
	int rho=0;
	sitio *SO = es->SO;
	sitio vecino;
	float Rho_P;
	int ON=es->ON;
	int NDX=es->NDX;
	int NDY=es->NDY;	
	int **s = es->s; 

	for(n=1;n<=ON;n++)
	{		
		for(w=1;w<=A;w++)
		{
		x=I_JKISS(0,diam) - radio;
		y=I_JKISS(0,diam) - radio;

		R2= x * x + y * y;
			if(R2<=r2 && R2!=0)
			{
				vecino.i=SO[n].i + x;
				vecino.j=SO[n].j + y;
				if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
				if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
				if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
				if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
			
				if( s[vecino.i][vecino.j]!=0 )
				{
					rho+=1;
				}
			}
		}
	}
	Rho_P=(float)rho/(float)ON;
return Rho_P; 	
}

float FuncionCorrelacion(estado *es,int radio)
{
	static float Rho_Ant=0.0;
	static int radio_Ant=0;
	float g,Rho;
	float Densidad=(float)(es->ON)/(float)((es->NDX)*(es->NDY));
	
	Rho=OnPromRadio(es, radio);
	
	if(radio_Ant != (radio-1))
	{
		radio_Ant=OnPromRadio(es, (radio-1));
	}
	
	g=(Rho-Rho_Ant)/(2.0*3.14159*radio*Densidad);
	
	Rho_Ant=Rho;
	radio_Ant=radio;
	
return g;
}

float FuncionCorrelacion2(estado *es,int radio)
{
	int x,y,n,i,j;
	int MasDer;
	sitio vecino;
	int **s=es->s;
	int Rho=0;
	int NDX=es->NDX;
	int NDY=es->NDY;
	int RadioInterior2 = (radio-1) * (radio-1);
	int Radio2 = radio * radio;
	int ON=es->ON;
	float Densidad=(float)ON/(float)(NDX*NDY);
	float g;
	
	for(n=1;n<=ON;n++)
	{
		x=es->SO[n].i;
		y=es->SO[n].j;
		MasDer=radio;
		i=MasDer;
		for(j=0;j<=i;j++)
		{ 
			do{
				if((i*i + j*j)<=Radio2)
				{
					vecino.i=x+i;				//Octante 1
					vecino.j=y+j;
					if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
					if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
					if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
					if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
						if(s[vecino.i][vecino.j]!=0){Rho++;}
						
					vecino.i=x-j;			//Octante 3
					vecino.j=y+i;
					if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
					if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
					if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
					if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
					if(s[vecino.i][vecino.j]!=0){Rho++;}	
							
					vecino.i=x-i;			//Octante 5
					vecino.j=y-j;
					if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
					if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
					if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
					if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
						if(s[vecino.i][vecino.j]!=0){Rho++;}	
						
					vecino.i=x+j;			//Octante 7
					vecino.j=y-i;
					if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
					if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
					if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
					if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
						if(s[vecino.i][vecino.j]!=0){Rho++;}
					
					if(i!=j && j!=0)
					{
						vecino.i=x+j;		//Octante 2
						vecino.j=y+i;
						if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
						if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
						if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
						if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
							if(s[vecino.i][vecino.j]!=0){Rho++;}

						vecino.i=x-i;			//Octante 4
						vecino.j=y+j;
						if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
						if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
						if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
						if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
							if(s[vecino.i][vecino.j]!=0){Rho++;}
												
									
						vecino.i=x-j;			//Octante 6
						vecino.j=y-i;
						if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
						if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
						if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
						if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
							if(s[vecino.i][vecino.j]!=0){Rho++;}
						
											
						vecino.i=x+i;			//Octante 8
						vecino.j=y-j;
						if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
						if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
						if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
						if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
							if(s[vecino.i][vecino.j]!=0){Rho++;}	
										
					}	
										
				}else{
					MasDer--;
				}
				i--;
			}while((i*i + j*j)>RadioInterior2);
			i=MasDer;
		}
	}
	
	g=(float)Rho/(6.283185307*radio*Densidad*ON);
	
return g;
}

float FuncionCorrelacionSpecies(estado *es,int radio,int TipoOrigen, int TipoDistante)
{
	int x,y,n,i,j;
	int MasDer;
	sitio vecino;
	int **s=es->s;
	int Rho=0;
	int NDX=es->NDX;
	int NDY=es->NDY;
	int RadioInterior2 = (radio-1) * (radio-1);
	int Radio2 = radio * radio;
	int ON=es->ON;
	float Densidad;
	float g;
	int **TIPO = es->TIPO;
	
	for(n=1;n<=ON;n++)
	{
		if(TIPO[es->SO[n].i][es->SO[n].j]==TipoOrigen)
		{
			x=es->SO[n].i;
			y=es->SO[n].j;
			MasDer=radio;
			i=MasDer;
			for(j=0;j<=i;j++)
			{ 
				do{
					if((i*i + j*j)<=Radio2)
					{
						vecino.i=x+i;				//Octante 1
						vecino.j=y+j;
						if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
						if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
						if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
						if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
							if(s[vecino.i][vecino.j]!=0 && TIPO[vecino.i][vecino.j]==TipoDistante){Rho++;}
							
						vecino.i=x-j;			//Octante 3
						vecino.j=y+i;
						if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
						if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
						if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
						if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
						if(s[vecino.i][vecino.j]!=0 && TIPO[vecino.i][vecino.j]==TipoDistante){Rho++;}	
								
						vecino.i=x-i;			//Octante 5
						vecino.j=y-j;
						if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
						if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
						if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
						if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
							if(s[vecino.i][vecino.j]!=0 && TIPO[vecino.i][vecino.j]==TipoDistante){Rho++;}	
							
						vecino.i=x+j;			//Octante 7
						vecino.j=y-i;
						if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
						if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
						if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
						if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
							if(s[vecino.i][vecino.j]!=0 && TIPO[vecino.i][vecino.j]==TipoDistante){Rho++;}
						
						if(i!=j && j!=0)
						{
							vecino.i=x+j;		//Octante 2
							vecino.j=y+i;
							if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
							if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
							if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
							if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
								if(s[vecino.i][vecino.j]!=0 && TIPO[vecino.i][vecino.j]==TipoDistante){Rho++;}

							vecino.i=x-i;			//Octante 4
							vecino.j=y+j;
							if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
							if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
							if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
							if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
								if(s[vecino.i][vecino.j]!=0 && TIPO[vecino.i][vecino.j]==TipoDistante){Rho++;}
													
										
							vecino.i=x-j;			//Octante 6
							vecino.j=y-i;
							if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
							if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
							if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
							if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
								if(s[vecino.i][vecino.j]!=0 && TIPO[vecino.i][vecino.j]==TipoDistante){Rho++;}
							
												
							vecino.i=x+i;			//Octante 8
							vecino.j=y-j;
							if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
							if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
							if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
							if(vecino.j > NDY){vecino.j = vecino.j - NDY;}	
								if(s[vecino.i][vecino.j]!=0 && TIPO[vecino.i][vecino.j]==TipoDistante){Rho++;}			
						}	
											
					}else{
						MasDer--;
					}
					i--;
				}while((i*i + j*j)>RadioInterior2);
				i=MasDer;
			}
		}
	}
	
	Densidad=(float)CuentaEspecie(es, TipoDistante)/(float)(NDX*NDY);
	g=(float)Rho/(6.283185307*((float)(radio*ON))*Densidad);
	
return g;
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
	
				if(s[i][j]==0){
				s[i][j]=1;
				TIPO[i][j]=tipo;
				(es->ON)++;
				SO[(es->ON)].i = i;
				SO[(es->ON)].j = j;
				INDICE[i][j]=(es->ON);
				}
return;
}

void ActualizaRhoVsT_MP(estado *es,Float2D_MP *RhoVsT,Dist_MP *Dist)	
{
int T=es->T;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;
int **TIPO=es->TIPO;
sitio *SO=es->SO;
float rho,rho_specie;
int NoEspecies = (RhoVsT->j_max);
int tot=NoEspecies + 1;
int rhoVec[tot];  
memset(rhoVec,0,tot * sizeof(int));

if(Dist!=NULL)
{
	float Dr=Dist->TamParticion;
}

int NoPart;

rho=((float)ON)/((float)(NDX*NDY));

#pragma omp atomic
RhoVsT->array[T][0]+=rho;	// Tipo 0 es el total

	if(NoEspecies!=0)
	{
	int n;
		for(n=1;n<=ON;n++)
		{
			rhoVec[TIPO[SO[n].i][SO[n].j]]+=1;
		}

		for(n=1;n<tot;n++)
		{
				if(rhoVec[n]!=0)
				{
					rho_specie=((float)rhoVec[n])/((float)(NDX*NDY));
					#pragma omp atomic
					RhoVsT->array[T][n]+=rho_specie;
					
					if(Dist!=NULL)
					{
						
						NoPart=(int)(rho_specie/Dr);
						#pragma omp atomic
						Dist->array[NoPart]++;
						Dist->T=T;
					}
				}
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
	int tam = (int)(1.0/Dist->TamParticion) + 1;
		Dist->array = (int *)calloc(tam, sizeof(int));
		if(imem==NULL)
		{
			puts("No se pudo alojar memoria para IntDist_MP\n");
			exit(0);
		}
		Dist->x_max=tam - 1;
		
return;
}

void JuntaFloat2D_MP(Float2D_MP *Origen, Float2D_MP *Destino)
{
int totj=(Destino->j_max);
int totx=Destino->x_max;
int T,tipo;
	for(T=0;T<=totx;T++)
	{
		for(tipo=0;tipo<=totj;tipo++)
		{	
		#pragma omp atomic
		Destino->array[T][tipo]+=Origen->array[T][tipo];
		}
	}
}

void JuntaDist_MP(Dist_MP *Origen, Dist_MP *Destino)
{
	
	int totx=Destino->x_max;
	int x;
		for(x=0;x<=totx;x++)
		{	
		#pragma omp atomic
		Destino->array[x]+=Origen->array[x];
		}

}


