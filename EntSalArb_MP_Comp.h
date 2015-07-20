/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */


void GuardaEstado(estado *es, FILE *archivo);
void CreaContenedor(char *nombre,runDescriptor run);
void GuardaEstadoEn(char *nombre, estado *es);
FILE* AbreRhoVsTEn(char *contenedor);
float ActualizaRhoVsT(estado *es,FILE *archivo,int NoEspecies);
int GuardaTiposEn(char *contenedor, estado *es);
FILE* AbreNoSpeciesVsTEn(char *contenedor);
void ActualizaNoSpeciesVsT(FILE *archivo,int Species, int T);
int CargaEstado(char *contenedor, char *nombre, estado *es, int NDX, int NDY);  //Usar sin haber alojado memoria antes!!!
void GuardaRhoVsT_MP(char *contenedor, Float2D_MP *RhoVsT, Dist_MP *RhoDist);
void GuardaTiposEn_MP(char *contenedor,Float2D_MP *MP_RhoVsT,int T);
void GuardaEstadoEn_MP(char *nombre, estado *es,int id,int ensamble);	
int CargaEstado_MP(char *contenedor, char *nombre, estado *es,int NDX, int NDY,int id,int NoEnsambles);  //Usar sin haber alojado memoria antes!!!
void GuardaCorrelacion_MP(char *contenedor, char *prefix, Float1D_MP *corr);
void GuardaCorrelacionTipo_MP(char *contenedor, Float1D_MP *corr);

void GuardaCorrXY(Float2D_MP *correlacion, char *contenedor, char *sufix);
void GuardaFloat1D_MP(char *contenedor,char *nombre, Float1D_MP *MP_Float1D);
void GuardaDist_MP(char *contenedor,char *nombre, Dist_MP *MP_Dist);
