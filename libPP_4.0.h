typedef struct {
int i;
int j;
} sitio;

typedef struct {
float Birth;
float Coagulation; 
float Dead;
int RadioBirth;
int RadioCoa;
} especie;

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

typedef struct {
float **array;
int i_max;
int j_max;
int NoEnsambles;
} Float2D_MP;

typedef struct {
int **array;
int i_max;
int j_max;
int NoEnsambles;
} Int2D_MP;

typedef struct {
int *array;
int T;
int x_max;
float TamParticion;
} Dist_MP;

typedef struct {
Float2D_MP RhoVsT;
Dist_MP Dist;
} RhoVsT_MP;

void SetBirth(float L, int tipo);
void SetCoagulation(float e, int tipo);
void SetDead(float d, int tipo);
void SetRadioBirth(int rb, int tipo);
void SetRadioCoa(int rc, int tipo);
void AlojaMemoria(estado *es, int NDX, int NDY);
void ResetEstado(estado *es);
void GeneraEstadoAleatorio(estado *es, float frac, int tipo);
void ActualizaRyC(estado *es, int N, int campo);
void BarrMCcRyCamp(estado *es);
void EligeUniforme(int i,int j,int radio, sitio *vecino);
void InsertaIndividuosAleatorio(estado *es, int N, int tipo);
float OnPromRadio(estado *es, int radio);
float FuncionCorrelacion(estado *es,int radio);
float FuncionCorrelacion2(estado *es,int radio);
float FuncionCorrelacionSpecies(estado *es,int radio,int TipoOrigen, int TipoDistante);
int CuentaEspecie(estado *es, int tipo);
void InsertaIndividuoEn(estado *es,int i,int j,int tipo);
void AlojaMemoriaEspecie(int tipo);
void EscalaTiempoMetabolico(int tipo);
void ActualizaRhoVsT_MP(estado *e,Float2D_MP *RhoVsT,Dist_MP *RhoDist);	
void IniciaMemoriaFloat2D_MP(Float2D_MP *ARRAY);
void IniciaMemoriaInt2D_MP(Int2D_MP *ARRAY);
void IniciaMemoriaDist_MP(Dist_MP *Dist);
void JuntaFloat2D_MP(Float2D_MP *Origen, Float2D_MP *Destino);
void JuntaDist_MP(Dist_MP *Origen, Dist_MP *Destino);
