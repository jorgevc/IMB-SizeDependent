#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "GNA.h"
#include "libPP_2.0.h"

void AlojaMemoria(estado *es, int NDX, int NDY)
{
int *imem;
int *INDmem;
sitio *SOmem;
int col, ile;
int k;
col = NDY + 1;
ile = NDX + 1;

	imem = (int *)malloc(col * ile * sizeof(int));
	INDmem = (int *)malloc(col * ile * sizeof(int));
	SOmem = (sitio *)malloc(col * ile * sizeof(sitio));
	if (imem == NULL || INDmem == NULL || SOmem == NULL)
	{	
		puts("\nNo se pudo alojar memoria para el estado\n");
		exit(0);
	}
	memset(imem,0,col * ile * sizeof(int));
	memset(INDmem,0,col * ile * sizeof(int));
	memset(SOmem,0,col * ile * sizeof(sitio));
	es->s = (int **)malloc(ile * sizeof(int *));
	es->INDICE = (int **)malloc(ile * sizeof(int *));
	es->SO = SOmem;
	if ((*es).s == NULL || (*es).INDICE == NULL)
	{	
		puts("\nNo se pudo alojar memoria para los apuntadores de la memoria\n");
		exit(0);
	}
	for (k = 0; k < ile; k++)
	{
		es->s[k] = imem + (k * col);
		es->INDICE[k] = INDmem + (k * col);
	}
	es->NDX = NDX;
	es->NDY = NDY;
	es->T = 0;
	es->ON = 0;
return;
}

void Actualiza4(estado *es, int N, float Lambda, float epsilon)
{
float Rand;
float pa, pcreacion;
int vecino;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i=es->SO[N].i;
int j=es->SO[N].j;
sitio *SO = es->SO;
	
	Rand = F_JKISS();
	//printf("Rand: %f\n",Rand); ///
	pa=1.0/(Lambda + 1.0 + epsilon);
	pcreacion = Lambda/(Lambda + 1.0 + epsilon);

	if(Rand<pa) //aniquilacion
	{
		s[i][j]=0;
		SO[N]=SO[(es->ON)];
		es->INDICE[SO[es->ON].i][SO[es->ON].j]=N;
		(es->ON)--;	
		
	}else{  //creation o coagulacion
		if(Rand<(pa + pcreacion)) //creation
		{
			vecino=I_JKISS(1,4);
		
			switch( vecino )
			{
				case 1:
					if(i != 1 && s[i-1][j]==0){ s[i-1][j]=1; (es->ON)++; 
						 SO[(es->ON)].i=(i-1);
						SO[(es->ON)].j=j;
						es->INDICE[(i-1)][j]=(es->ON);
					}
					else{
					//	printf("NDX= %d, s[NDX][j] = %d\n",NDX,s[NDX][j]);    //borrar
						if(i == 1 && s[NDX][j]==0){s[NDX][j]=1; (es->ON)++;
							 SO[(es->ON)].i=NDX;
							SO[(es->ON)].j=j;
							es->INDICE[NDX][j]=(es->ON);
						 }
					}
				break;
				case 2:
					if(i != NDX && s[i+1][j]==0){ s[i+1][j]=1; (es->ON)++;
						 SO[(es->ON)].i=(i+1);
						SO[(es->ON)].j=j;
						es->INDICE[(i+1)][j]=(es->ON);
					 }
					else{
						if(i == NDX && s[1][j]==0){s[1][j]=1; (es->ON)++; 
							 SO[(es->ON)].i=1;
							SO[(es->ON)].j=j;
							es->INDICE[1][j]=(es->ON);
						}
					}
				break;
				case 3:
					if(j != 1 && s[i][j-1]==0){ s[i][j-1]=1; (es->ON)++;
						 SO[(es->ON)].i=i;
						SO[(es->ON)].j=(j-1);
						es->INDICE[i][j-1]=(es->ON);
					 }
					else{
						if(j == 1 && s[i][NDY]==0){s[i][NDY]=1; (es->ON)++;
							 SO[(es->ON)].i=i;
							SO[(es->ON)].j=NDY;
							es->INDICE[i][NDY]=(es->ON);
						 }
					}
				break;
				case 4:
					if(j != NDY && s[i][j+1]==0){ s[i][j+1]=1; (es->ON)++;
						 SO[(es->ON)].i=i;
						SO[(es->ON)].j=(j+1);
						es->INDICE[i][j+1]=(es->ON);
					 }
					else{
						if(j == NDY && s[i][1]==0){s[i][1]=1; (es->ON)++;
							 SO[(es->ON)].i=i;
							SO[(es->ON)].j=1;
							es->INDICE[i][1]=(es->ON);
						 }
					}
				break;
			}
		}else{ //coagulacion
			vecino=I_JKISS(1,4);
		
			switch( vecino )
			{
				case 1:
					if(i != 1 && s[i-1][j]==1){ s[i-1][j]=0;  
						SO[(es->INDICE[i-1][j])]=SO[(es->ON)];
						es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[i-1][j]);
						(es->ON)--;
					}
					else{
						if(i == 1 && s[NDX][j]==1){s[NDX][j]=0; 
							SO[(es->INDICE[NDX][j])]=SO[(es->ON)];
							es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[NDX][j]);
							(es->ON)--;
						 }
					}
				break;
				case 2:
					if(i != NDX && s[i+1][j]==1){ s[i+1][j]=0;
						SO[(es->INDICE[i+1][j])]=SO[(es->ON)];
						es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[i+1][j]);
						 (es->ON)--;
					 }
					else{
						if(i == NDX && s[1][j]==1){s[1][j]=0;
							SO[(es->INDICE[1][j])]=SO[(es->ON)];
							es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[1][j]);
							 (es->ON)--; 
						}
					}
				break;
				case 3:
					if(j != 1 && s[i][j-1]==1){ s[i][j-1]=0;
						SO[(es->INDICE[i][j-1])]=SO[(es->ON)];
						es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[i][j-1]);
						 (es->ON)--;
					 }
					else{
						if(j == 1 && s[i][NDY]==1){s[i][NDY]=0; 
							SO[(es->INDICE[i][NDY])]=SO[(es->ON)];
							es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[i][NDY]);
							(es->ON)--;
						 }
					}
				break;
				case 4:
					if(j != NDY && s[i][j+1]==1){ s[i][j+1]=0;
						SO[(es->INDICE[i][j+1])]=SO[(es->ON)];
						es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[i][j+1]);
						 (es->ON)--;
					 }
					else{
						if(j == NDY && s[i][1]==1){s[i][1]=0;
							SO[(es->INDICE[i][1])]=SO[(es->ON)];
							es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[i][1]);
							 (es->ON)--;
						 }
					}
				break;
			}

		}
	}
return;
}

void BarridoMC4(estado *es, float Lambda, float Epsilon)
{
int Indice;
float DT=0.0;

	while(DT<1.0){
		if((es->ON) != 0){
			 DT+=1.0/(es->ON); 
			Indice = I_JKISS(1,(es->ON));
			Actualiza4(es,Indice,Lambda,Epsilon);	
		}else{
			DT=2.0;
		}
	}

(es->T)++;				
		
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
RellenaIndiceYSO(es);
es->T=0;
return;
}

void RellenaIndiceYSO(estado *es)
{
int i,j;
int NDX = es->NDX;
int NDY = es->NDY;
int **s = es->s;
int N=0;
sitio *SO = es->SO;
int **INDICE = es->INDICE;

	for(i=1;i<=NDX;i++)
	{
		for(j=1;j<=NDY;j++)
		{
			if(s[i][j]==1)
			{
				N++;
				SO[N].i = i;
				SO[N].j = j;
				INDICE[i][j]=N;
			}
		}
	}
return;
}

void ActualizaR(estado *es, int N, float Lambda, float epsilon, float dead, int radioCre, int radioCua)
{
float Rand;
float pa, pcreacion, C;
sitio vecino;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i=es->SO[N].i;
int j=es->SO[N].j;
sitio *SO = es->SO;
	
	Rand = F_JKISS();
	//printf("Rand: %f\n",Rand); ///
	
	C=dead+Lambda+epsilon;
	pa=dead/C;
	pcreacion = Lambda/C;

	if(Rand<pa) //aniquilacion
	{
		s[i][j]=0;
		SO[N]=SO[(es->ON)];
		es->INDICE[SO[es->ON].i][SO[es->ON].j]=N;
		(es->ON)--;	
		
	}else{  //creation o coagulacion
		if(Rand<(pa + pcreacion)) //creation
		{
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
			}
			
		}else{ //coagulacion
			
			 EligeUniforme(i,j,radioCua,&vecino);
			
			if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
			if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
			if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
			if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
			
			if(s[vecino.i][vecino.j]==1)
			{
				s[vecino.i][vecino.j]=0;
				SO[(es->INDICE[vecino.i][vecino.j])]=SO[(es->ON)];
				es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[vecino.i][vecino.j]);
				(es->ON)--;
			}
		}
	}
return;
}

void BarridoMCR(estado *es, float Lambda, float Epsilon, float dead, int radioCre, int radioCua)
{
int Indice;
float DT=0.0;


	while(DT<1.0){
		if((es->ON) != 0){
			 DT+=1.0/(es->ON); 
			Indice = I_JKISS(1,(es->ON));
			ActualizaR(es,Indice,Lambda,Epsilon,dead,radioCre,radioCua);	
		}else{
			DT=2.0;
		}
	}

(es->T)++;				
		
return;
}

void EligeUniforme(int i,int j,int radio, sitio *vecino)
{
	int x,y,R2;
	int diam = 2 * radio;
	int r2 = radio * radio;
	
		do{
		x=I_JKISS(0,diam) - radio;
		y=I_JKISS(0,diam) - radio;

		R2= x * x + y * y;
		}while( R2 > r2 || R2==0)
		
		vecino->i = i + x;
		vecino->j = j + y;

return; 	
}

