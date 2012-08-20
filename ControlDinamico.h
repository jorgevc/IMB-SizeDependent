
#include "InterfaseControlDinamico.h"

int* activa_escucha(int clave);
void cierra_escucha(void);
void Servicio(int T,char *contenedor,Float2D_MP *MP_RhoVsT, Float2D_MP *MP_RhoVsT_1, Dist_MP *MP_RhoDist, Dist_MP *MP_RhoDist_1, estado *e, int MaxPar);
void SalidaCD(int *i,int T_max);
