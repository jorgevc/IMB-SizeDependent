#include <stdio.h>
#include <string.h>
#include "GNA.h"


main(){
FILE *dist;
int i,dado;

dist=fopen("borrar","w");

	for(i=1;i<=1000;i++)
	{
		dado=I_JKISS(1,1000);
		fprintf(dist,"%d %d\n",i,dado);
	}
	
	fclose(dist);
}
