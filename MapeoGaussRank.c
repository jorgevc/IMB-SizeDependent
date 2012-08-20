#include <stdio.h>
#include <math.h>

main()
{
	FILE *ARCH;
	
double xIni = 0.0; //exp(-4.64);
double xFin = 0.02; //exp(-4.59);
double Dx = (xFin - xIni)/400.0;

double C = 0.0;
double xo=exp(-4.61); //0.00982556;  //exp(-4.61);
double s = 0.000104243; // (xIni-xFin) 0.001734243;  0.000134243
double x;

printf("%1.10f, %1.8lf",s,s);
	for(x=xFin;x>=xIni;x-=Dx)
	{
		C+=exp(-(pow(x-xo,2.0))/(2*pow(s,2.0)));
	}
C = 50.0/C;


 ARCH = fopen("DATOS/MapeoGaussRank.dat","w");

double y1,g;
double tmp;
int x1 = 1;
for(x=xFin;x>=xIni;x-=Dx)
{
	g=C*exp(-(pow(x-xo,2.0))/(2.0*pow(s,2.0)));
	y1=x;
	tmp=g;
	while(tmp>=0.5)
	{
		fprintf(ARCH,"%d %lf %lf %lf\n",x1,y1,x,g);
		x1++;
		y1 -= Dx/g;
		tmp--;
	}
}

fclose(ARCH);

return;
}
