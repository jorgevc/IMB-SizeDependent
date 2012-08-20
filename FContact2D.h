typedef struct {
int NDX;
int NDY;
int T;
int **s;
int ON;
} estado;

typedef struct {
int i;
int j;
} sitio;

void AlojaMemoria(estado * es, int NDX, int NDY);
int CuentaVecinos(estado * es, int i, int j);
void Actualiza(estado * es, int i,int j, float Lambda);
void BarridoMonteCarlo(estado * es, float Lambda);
void GeneraEstadoAleatorio(estado * es, float frac);
void Actualiza2(estado *es,int i, int j, float Lambda);
void BarridoMC2(estado *es,float Lambda);
void BarridoMC3(estado *es, sitio *SO, float Lambda);
void Actualiza3(estado *es,sitio *SO,int N, float Lambda);
void RellenaSO(estado *es,sitio *SO);
