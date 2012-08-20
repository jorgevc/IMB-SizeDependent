#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
//#include "FContact2D.h"
//#include "libPP_2.0.h"
#include "libPP_3.0.h"
#include "EntSal2D.h"

void GuardaEstado(estado *es, FILE *archivo)
{
int **s=es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i,j;


	for(i=1;i<=NDX;i++)
	{
		for(j=1;j<=NDY;j++)
		{
			if(s[i][j]==1){ fprintf(archivo,"%d   %d\n",i,j); }
			//else{fprintf(archivo,"No hay datos para: %d   %d \n",i,j);}
		}
	}
}

FILE* OpenFile(char *nombre,int T,float rho)
{
FILE *aA;
char archivo[100]="DATOS/";
// mkdir(archivo,00777); Descomentar para crear el directorio si es que no existe.
strcat(archivo,nombre);
		aA=fopen(archivo, "w");
		fprintf(aA,"\n# x   y   rho=%f   T=%d   \n",rho,T);
		return aA;
}

void CreaContenedor(char *nombre)
{
char dir[100]="DATOS/";

strcat(dir,nombre);
mkdir(dir,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH));
return;
}

void GuardaEstadoEn(char *nombre, estado *es)
{
int T=es->T;
char paso[15];
char archivo[100];
FILE *datos;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;
float rho;

rho=((float)ON)/((float)(NDX*NDY));


sprintf(paso,"/T_%03d",T);
strcpy(archivo,nombre);
strcat(archivo,paso);

	datos=OpenFile(archivo,T,rho);
	GuardaEstado(es, datos);
	fclose(datos);

return;
}

FILE* AbreRhoVsTEn(char *nombre)
{
FILE *aA;
char archivo[100]="DATOS/";
// mkdir(archivo,00777); Descomentar para crear el directorio si es que no existe.
strcat(archivo,nombre);
strcat(archivo,"/RhoVsT");
		aA=fopen(archivo, "w");
		fputs("# t   rho \n",aA);
		fclose(aA);
		aA=fopen(archivo, "a");
		return aA;
}

float ActualizaRhoVsT(estado *es,FILE *archivo)
{
int T=es->T;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;
float rho;

rho=((float)ON)/((float)(NDX*NDY));

fprintf(archivo,"%d   %f\n",T,rho);
return rho;
}
