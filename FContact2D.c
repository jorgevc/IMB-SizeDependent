/*
Copyright 2012 Jorge Velazquez
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "FContact2D.h" 
#include "GNA.h"

void AlojaMemoria(estado *es, int NDX, int NDY)
{
int *imem;
int col, ile;
int k;
col = NDY + 1;
ile = NDX + 1;

	imem = (int *)malloc(col * ile * sizeof(int));
	if (imem == NULL)
	{	
		puts("\nNo se pudo alojar memoria para el estado\n");
		exit(0);
	}
	memset(imem,0,col * ile * sizeof(int));
	es->s = (int **)malloc(ile * sizeof(int *));
	if ((*es).s == NULL)
	{	
		puts("\nNo se pudo alojar memoria para los apuntadores de la memoria\n");
		exit(0);
	}
	for (k = 0; k < ile; k++)
	{
		es->s[k] = imem + (k * col);
	}
	es->NDX = NDX;
	es->NDY = NDY;
	es->T = 0;
	es->ON = 0;
	
return;
}

int CuentaVecinos(estado *es, int i, int j)
{
int NDX=es->NDX;
int NDY=es->NDY;
int **s;
s=es->s;

int n=0; // numero de vecionos
	if(i != 1 && i != NDX && j != 1 && j != NDY)
	{
		if(s[i-1][j]==1){n++;}
		if(s[i+1][j]==1){n++;}
		if(s[i][j-1]==1){n++;}
		if(s[i][j+1]==1){n++;}
	}else{
		if(i==1 || i==NDX){
			if(i==1){
				if(s[NDX][j]==1){n++;}
				if(s[2][j]==1){n++;}
			}
			if(i==NDX){
				if(s[i-1][j]==1){n++;}
				if(s[1][j]==1){n++;}
			}
		}else{
			if(s[i-1][j]==1){n++;}
			if(s[i+1][j]==1){n++;}
		}
		if(j==1 || j==NDY){		
			if(j==1){
				if(s[i][NDY]==1){n++;}
				if(s[i][2]==1){n++;}
			}
			if(j==NDY){
				if(s[i][j-1]==1){n++;}
				if(s[i][1]==1){n++;}
			}
		}else{
			if(s[i][j-1]==1){n++;}
			if(s[i][j+1]==1){n++;}
		}
	}

return n;
}

void Actualiza(estado *es, int i, int j, float Lambda)
{	
	float dimension = 2.0;
	float p;
	float vecinos;
	int **s;
	s=es->s;
	float Rand;

	if(s[i][j]==0)
	{
		Rand = F_JKISS();
	
	 	vecinos=(float)CuentaVecinos(es,i,j);
	
		p=Lambda*vecinos/((2.0)*dimension);

		if (Rand<p)
		{
			s[i][j]=1;
			(es->ON)++;
		}
	}else{
		if(s[i][j]==1)
		{
			s[i][j]=0;
			(es->ON)--;
		}
	}
return;
}

void BarridoMonteCarlo(estado * es, float Lambda)
{
int NDX;
int NDY;
NDX=es->NDX;
NDY=es->NDY;
int N;
int i;
N=NDX*NDY;
int RandX;
int RandY;
long int xrand;
long int yrand;
	for (i=1;i<=N;i++){
		RandX = I_JKISS(1, NDX);
		RandY = I_JKISS(1, NDY);
		Actualiza(es,RandX,RandY,Lambda);
	}
es->T++;
return;
}

void GeneraEstadoAleatorio(estado *es, float frac)
{
int NDX = es->NDX;
int NDY = es->NDY;
float fNDX = es->NDX;
float fNDY = es->NDY;
int col = es->NDY + 1;
int ile = es->NDX + 1;
int **s = es->s; 
long int xrand;
long int yrand;
int n=0;
int N;
int i,j;

	if(frac==0)
	{
		memset(*s,0,col * ile * sizeof(int));
		(es->ON)=0;
	}
	if(frac==1)
	{
		(es->ON)=0;
		for(i=1;i<=NDX;i++)
		{
			for(j=1;j<=NDY;j++)
			{
			s[i][j]=1;
			(es->ON)++;
			}
		}
	}
	if(frac <= 0.5 && frac != 0.0)
	{
		N = (int)((fNDX * fNDY) * frac);
		memset(*s,0,col * ile * sizeof(int));
		while(n<N)
		{ 
			i = I_JKISS(1, NDX);
			j = I_JKISS(1, NDY);		
			if(s[i][j]==0){
				s[i][j]=1;
				n++;
			}
		}
		(es->ON)=N;	
	}
	if(frac > 0.5 && frac != 1.0)
	{
		N = (int)((fNDX * fNDY) * (1.0-frac));
		(es->ON)=0;
		for(i=1;i<=NDX;i++)
		{
			for(j=1;j<=NDY;j++)
			{
			s[i][j]=1;
			(es->ON)++;
			}
		}
		while(n<N)
		{
			i = I_JKISS(1, NDX);
			j = I_JKISS(1, NDY);			
			if(s[i][j]==1){
				s[i][j]=0;
				n++;
				(es->ON)--;
			}
		}
	}
es->T=0;
return;
}

void PasoMC2(estado *es, float Lambda)
{
int Indice,i,j;
int CuentaON = 0;
int ON = es->ON;
int NDY = es->NDY;
int NDX = es->NDX;
int **s = es->s;

	Indice = I_JKISS(1,ON);

	// printf("NDX=%d ,NDY=%d\n",NDX,NDY);   //borrar

	for (i=1;i<=NDX;i++)
	{
		for(j=1;j<=NDY;j++)
		{
			if(s[i][j]==1)
			{
				CuentaON++;
				if(CuentaON == Indice){
				//	printf("i=%d, j=%d\n",i,j);   //borrar
					Actualiza2(es,i,j,Lambda);
					return;
				}
			}
		}
	}
return;
}

void Actualiza2(estado *es,int i, int j, float Lambda)
{
float Rand;
float pa;
int vecino;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
	
	Rand = F_JKISS();

	pa=1.0/(Lambda + 1.0);

	if(Rand<pa) //aniquilacion
	{
		s[i][j]=0;
		(es->ON)--;			
	}else{  //creation
		vecino=I_JKISS(1,4);
	//	printf("vecino: %d , valor i=%d, valor j=%d\n",vecino,i,j);   //borrrar
		switch( vecino )
		{
			case 1:
				if(i != 1 && s[i-1][j]==0){ s[i-1][j]=1; (es->ON)++; }
				else{
				//	printf("NDX= %d, s[NDX][j] = %d\n",NDX,s[NDX][j]);    //borrar
					if(i == 1 && s[NDX][j]==0){s[NDX][j]=1; (es->ON)++; }
				}
			break;
			case 2:
				if(i != NDX && s[i+1][j]==0){ s[i+1][j]=1; (es->ON)++; }
				else{
					if(i == NDX && s[1][j]==0){s[1][j]=1; (es->ON)++; }
				}
			break;
			case 3:
				if(j != 1 && s[i][j-1]==0){ s[i][j-1]=1; (es->ON)++; }
				else{
					if(j == 1 && s[i][NDY]==0){s[i][NDY]=1; (es->ON)++; }
				}
			break;
			case 4:
				if(j != NDY && s[i][j+1]==0){ s[i][j+1]=1; (es->ON)++; }
				else{
					if(j == NDY && s[i][1]==0){s[i][1]=1; (es->ON)++; }
				}
			break;
		}
	}
return;
}

void BarridoMC2(estado *es,float Lambda)
{
int ON = (es->ON);  //analizar mejor ??
// printf("Iteraciones en el paso:%d\n",ON);
int i;  
	for (i=1;i<=ON;i++){
		PasoMC2(es,Lambda);
	}
es->T++;
return;
}

void Actualiza3(estado *es,sitio *SO,int N, float Lambda)
{
float Rand;
float pa;
int vecino;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i=SO[N].i;
int j=SO[N].j;

//if(Lambda>1.648){printf("Indice: %d, ON: %d\n",N,(es->ON));} //borrar

	
	Rand = F_JKISS();

	pa=1.0/(Lambda + 1.0);

	if(Rand<pa) //aniquilacion
	{
//if(Lambda>1.648){printf("Se aniquila i=%d , j=%d \n",i,j);} //borrar
		s[i][j]=0;
		SO[N]=SO[(es->ON)];
		(es->ON)--;			
	}else{  //creation
		vecino=I_JKISS(1,4);
	//	printf("vecino: %d , valor i=%d, valor j=%d\n",vecino,i,j);   //borrrar
	//if(Lambda>1.648){  printf("Se contagia i=%d , j=%d \n",i,j);} //borrar
		switch( vecino )
		{
			case 1:
				if(i != 1 && s[i-1][j]==0){ s[i-1][j]=1; (es->ON)++; 
					 SO[(es->ON)].i=(i-1);
					SO[(es->ON)].j=j;
				}
				else{
				//	printf("NDX= %d, s[NDX][j] = %d\n",NDX,s[NDX][j]);    //borrar
					if(i == 1 && s[NDX][j]==0){s[NDX][j]=1; (es->ON)++;
						 SO[(es->ON)].i=NDX;
						SO[(es->ON)].j=j;
					 }
				}
			break;
			case 2:
				if(i != NDX && s[i+1][j]==0){ s[i+1][j]=1; (es->ON)++;
					 SO[(es->ON)].i=(i+1);
					SO[(es->ON)].j=j;
				 }
				else{
					if(i == NDX && s[1][j]==0){s[1][j]=1; (es->ON)++; 
						 SO[(es->ON)].i=1;
						SO[(es->ON)].j=j;
					}
				}
			break;
			case 3:
				if(j != 1 && s[i][j-1]==0){ s[i][j-1]=1; (es->ON)++;
					 SO[(es->ON)].i=i;
					SO[(es->ON)].j=(j-1);
				 }
				else{
					if(j == 1 && s[i][NDY]==0){s[i][NDY]=1; (es->ON)++;
						 SO[(es->ON)].i=i;
						SO[(es->ON)].j=NDY;
					 }
				}
			break;
			case 4:
				if(j != NDY && s[i][j+1]==0){ s[i][j+1]=1; (es->ON)++;
					 SO[(es->ON)].i=i;
					SO[(es->ON)].j=(j+1);
				 }
				else{
					if(j == NDY && s[i][1]==0){s[i][1]=1; (es->ON)++;
						 SO[(es->ON)].i=i;
						SO[(es->ON)].j=1;
					 }
				}
			break;
		}
	}
return;
}

void BarridoMC3(estado *es, sitio *SO, float Lambda)
{
int Indice;
float DT=0.0;

	while(DT<1.0){
		if((es->ON) != 0){
			 DT+=1.0/(es->ON); 
			Indice = I_JKISS(1,(es->ON));
			Actualiza3(es,SO,Indice,Lambda);	
		}else{
			DT=2.0;
		}
	}

(es->T)++;				
		
return;
}

void RellenaSO(estado *es,sitio *SO)
{
int i,j;
int NDX = es->NDX;
int NDY = es->NDY;
int **s = es->s;
int N=0;

	for(i=1;i<=NDX;i++)
	{
		for(j=1;j<=NDY;j++)
		{
			if(s[i][j]==1)
			{
				N++;
				SO[N].i = i;
				SO[N].j = j;
			}
		}
	}
return;
}


