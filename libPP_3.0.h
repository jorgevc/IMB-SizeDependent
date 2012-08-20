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
int **TIPO;
} estado;

void SetBirth(float L);
void SetCoagulation(float e);
void SetDead(float d);
void SetRadioBirth(int rb);
void SetRadioCoa(int rc);
void AlojaMemoria(estado *es, int NDX, int NDY);
void ResetEstado(estado *es);
void GeneraEstadoAleatorio(estado *es, float frac, int tipo);
void ActualizaRyC(estado *es, int N, int campo);
void BarrMCcRyCamp(estado *es);
//float campo(int i, int j)
void EligeUniforme(int i,int j,int radio, sitio *vecino);
void InsertaIndividuosAleatorio(estado *es, int N, int tipo);
float OnPromRadio(estado *es, int radio);
float FuncionCorrelacion(estado *es,int radio);
float FuncionCorrelacion2(estado *es,int radio);
float FuncionCorrelacionSpecies(estado *es,int radio,int TipoOrigen, int TipoDistante);
int CuentaEspecie(estado *es, int tipo);
void InsertaIndividuoEn(estado *es,int i,int j,int tipo);
