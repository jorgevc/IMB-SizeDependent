#define MAXNUCLEOS 15

typedef struct{
int servidor_ocupado;
int guarda_estado;
int guarda_RhoVsT;
int mostrar_tiempos;
int salida_loop;
int salida_violenta;
int tiempos[MAXNUCLEOS][2];
char mensaje[100];
} interfase;
