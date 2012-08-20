void GuardaEstado(estado *es, FILE *archivo);
FILE* OpenFile(char *nombre, int T, float rho);
void CreaContenedor(char *nombre);
void GuardaEstadoEn(char *nombre, estado *es);
FILE* AbreRhoVsTEn(char *nombre);
float ActualizaRhoVsT(estado *es,FILE *archivo);
