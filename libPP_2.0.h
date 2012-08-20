typedef struct {
int i;
int j;
} sitio;

typedef struct {
int NDX;
int NDY;
int T;
int ON;
int **s; 
int **INDICE;
sitio *SO;
} estado;

void AlojaMemoria(estado * es, int NDX, int NDY);
void Actualiza4(estado *es, int N, float Lambda, float epsilon);
void BarridoMC4(estado *es, float Lambda, float Epsilon);
void GeneraEstadoAleatorio(estado * es, float frac);
void RellenaIndiceYSO(estado *es);
void ActualizaR(estado *es, int N, float Lambda, float epsilon, float dead, int radioCre, int radioCua);
void EligeUniforme(int i,int j,int radio, sitio *vecino);
