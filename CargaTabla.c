#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
//#include "libPP_5.0.h"
//#include "EntSalArb_MP.h"


main()
{
	
float TamParticion=2;
float xFin=250.0;
float xIni=0.0;

	float Dx = xFin - xIni;
	int i_max = (int)(Dx/TamParticion);
	
float tmp[i_max + 1][2];


//int i_max=500;
//float corr1[501];
//float corr2[501];

int i;
for(i=0;i<=i_max;i++)
{
	//corr1[i]=0.0;
	//corr2[i]=0.0;
	tmp[i][0]=0.0;
	tmp[i][1]=0.0;
}
	
FILE *datos=NULL;
FILE *Arch=NULL;
size_t tam_buffer=100*sizeof(char);
char *buffer;
buffer = (char *) malloc (tam_buffer + 1);
int args_assigned = 0;
char contenedor[200] = "Crecimiento2_(B,D,C,RB,RC)@(0.000,0.000,2.0000,10,60)_(NDX,Tmax)@(500,10000)";
char nombre[50] = "Mark_CorrT_2201";
char nombre2[50] = "";
char nombreRecipiente[50] = "Mark_CorrT_2201.smt";

char recipiente[250];
char archivo[250];
char archivo2[250];

sprintf(archivo,"DATOS/%s/%s",contenedor,nombre);
sprintf(archivo2,"DATOS/%s/%s",contenedor,nombre2);
sprintf(recipiente,"DATOS/%s/%s",contenedor,nombreRecipiente);

float x;
float pr;
int Nuevo_i;

	if((datos = fopen (archivo, "r"))==NULL){
		puts("\nNo se pudo abrir para leer1\n");
		return 0;
		}
	
	while(getline(&buffer, &tam_buffer, datos)!=-1)
	{
		if(strchr(buffer, '#')==NULL)
		{
			args_assigned = sscanf(buffer, "%*s %s %f %*d %f %f", &spec, &size, &x , &y);
			if(args_assigned==4)
			{
				//corr1[x]=pr;
				
				Nuevo_i=(int)((x-xIni)/TamParticion);
				Nuevo_j=(int)((y-yIni)/TamParticion);
				if(Nuevo_i>=0 && Nuevo_i<=i_max && Nuevo_j>=0 && Nuevo_j<=j_max)
				{
					tmp[Nuevo_i][Nuevo_j]++;
					//tmp[Nuevo_i][1]=tmp[Nuevo_i][1] + 1.0;
				//	printf("leido %d valor:%f\n",Nuevo_i,pr);
				}else{
					printf("No cupo i=%d\n",Nuevo_i);
				}			 
			}
		}
	}
	
	fclose(datos);
	/////////////////
	//if((datos = fopen (archivo2, "r"))==NULL){
		//puts("\nNo se pudo abrir para leer2\n");
		//return 0;
		//}
	
	//while(getline(&buffer, &tam_buffer, datos)!=-1)
	//{
		//if(strchr(buffer, '#')==NULL)
		//{
			//args_assigned = sscanf(buffer, "%d %f", &x, &pr);
			//if(args_assigned==2 && x<=500)
			//{
				//corr2[x]=pr;
			//}
		//}
	//}
	
	//fclose(datos);
	//////////////
	//float Integral = 0.0;
	//for(i=1;i<=250;i++)
	//{
		//Integral += (corr1[i]-corr2[i]);
	//}
	
	//printf("Resutlado de %s \nNueva Medida = %f\n",contenedor,Integral);
	
	/////////////////////////////////////
	fftw_complex *out;
	fftw_complex *in;
	fftw_plan p;
	int NDX = i_max;
	int NDY = j_max;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * NDY);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NDX * NDY);

	p = fftw_plan_dft_2d(NDX * NDY, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	 for ( i = 0; i < NDX; i++ )
	  {
		  in[i][0] = (corr1[i+1]);
		  in[i][1] = 0.0;
	  }
		
	//fftw_execute(p); 
	
	
	
	if((Arch = fopen (recipiente, "w"))==NULL){
		puts("\nNo se pudo abrir para leer3\n");
		return 0;
		}
		
		//for(i=0;i<i_max;i++)
		//{
			//fprintf(Arch,"%d %f %f\n",i,out[i][0]/corr1[1] + 1.0, out[i][1]/corr1[1] + 1.0);
		//}
	////////////////////////////////////
	for(i=0;i<=i_max;i++)
	{
		//if(tmp[i][1]!=0.0)
		//{
				fprintf(Arch,"%f %f\n",(xIni + ((float)i) * (TamParticion)) ,tmp[i][0]/tmp[i][1]);
		//}
	}
	////////////////////////////////////
	fclose(Arch);
	
	
	
return;
}
