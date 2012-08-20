#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "FContact2D.h" 
#include "GNA.h"




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

