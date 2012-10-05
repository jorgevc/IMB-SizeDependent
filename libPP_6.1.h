typedef struct {
int i;
int j;
} sitio;

typedef struct {
float Birth;
float Coagulation; 
float Dead;
float CoagulationIntra;
int RadioBirth;
int RadioCoa;
int RadioCoaIntra;
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
float Meta_T;
float Max_Metabolic;
int **AGE;
} estado;

typedef struct {
float **array;
int i_max;
int j_max;
int T;
int NoEnsambles;
} Float2D_MP;

typedef struct {
float *array;
int i_max;
int T;
int NoEnsambles;
} Float1D_MP;

typedef struct {
int **array;
int i_max;
int j_max;
int NoEnsambles;
} Int2D_MP;

typedef struct {
int *array;
int T;
int i_max;
float TamParticion;
int NoEnsambles;
float xIni;
float xFin;
} Dist_MP;

typedef struct {
float CoaFact;
float CoaExp;
float MetFact;
float MetExp;
float ResurcesFact;	
} model;

void SetBirth(float L, int tipo);
void SetCoagulation(float e, int tipo);
void SetCoagulationIntra(float e,int tipo);
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
float FuncionCorrelacion(estado *es,int radio);			//No eficiente
float FuncionCorrelacion2(estado *es,int radio);		//No eficiente
float FuncionCorrelacionSpecies(estado *es,int radio,int TipoOrigen, int TipoDistante); //No eficiente
int CuentaEspecie(estado *es, int tipo);
void InsertaIndividuoEn(estado *es,int i,int j,int tipo);
void AlojaMemoriaEspecie(int tipo);
void EscalaTiempoMetabolico(int tipo);
void ActualizaRhoVsT_MP(estado *e,Float2D_MP *RhoVsT,Dist_MP *RhoDist);	
void IniciaMemoriaFloat2D_MP(Float2D_MP *ARRAY);
void IniciaMemoriaInt2D_MP(Int2D_MP *ARRAY);
void IniciaMemoriaDist_MP(Dist_MP *Dist);
void SumaFloat2D_MP(Float2D_MP *Origen, Float2D_MP *Destino);
void SumaDist_MP(Dist_MP *Origen, Dist_MP *Destino);
void InicializaFloat2D_MP(Float2D_MP *Objeto, int i_max, int j_max, int NoEnsambles);
void InicializaDist_MP(Dist_MP *Objeto, float TamParticion);
void SetSpecie(int NoEspecie, float Birth, float Coagulation, float Dead, float RadioBirth, float RadioCoa);
void ResetFloat2D_MP(Float2D_MP *ARRAY);
void ResetDist_MP(Dist_MP *Dist);
void SetSpecie2(int NoEspecie, float Birth, float Coagulation, float CoagulationIntra, float Dead, float RadioBirth, float RadioCoa, float RadioCoaIntra);
void ActualizaCorrelacion_MP(estado *es, Float1D_MP *corr);
void ActualizaCorrelacionTipo_MP(estado *es, Float1D_MP *corr, int TipoOrigen, int TipoObjetivo);
void SumaFloat1D_MP(Float1D_MP *Origen,Float1D_MP *Destino);
void InicializaFloat1D_MP(Float1D_MP *Objeto, int i_max);

void GeneraEstadoAleatorioTamano(estado *es, float frac, int tipo, int tamano);
void ActualizaRyCTamano(estado *es, int N, int campo);
void BarrMCcRyCampTamano(estado *es, float flujo_recursos, model *param);
//void InsertaIndividuosAleatorioTamano(estado *es, int N, int tipo, int tamano);
//float CorrelacionEspeciesTamano(estado *es,int radio,int TamOrigen, int TamDistante);
void ActualizaDistTamano_MP(estado *e, Float1D_MP *TamDist, char Opcion);
//void ActualizaCorrelacionTamano_MP(estado *es, Float1D_MP *corr, int TamanoOrigen, int TamanoObjetivo);
void ActualizaRecursos_MP(estado *es,Float2D_MP *RhoVsT);

void CFFT_Tipos_MP(estado *es, int NoEnsambles, Float2D_MP *correlacion,int TipoOrigen, int TipoDestino);
void CFFT_Mark_MP(estado *es, int NoEnsambles, Float2D_MP *correlacion,int TamanoOrigen, int TamanoDestino);
void DoblaCorrelacion(Float2D_MP *corr2D);
void CompactaCorrelacion(Float2D_MP *corr2D, Float1D_MP *corrRadial);
void ResetFloat1D_MP(Float1D_MP *ARRAY);

void ActualizaUniv(estado *es, int N, model *modelo);
void LiberaMemoria(estado *es);

/** Calculates the likelyhood of Origin with Experiment 
 * 
 */
float LikelyHood(Float1D_MP *Origin, Float1D_MP *Experiment);
void CargaExperiment(Float1D_MP *Experiment);
