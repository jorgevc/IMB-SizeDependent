PDSpaceStructure: PDSpaceStructure.c libPP_5.0.c EntSalArb_MP.c GNA.c ControlDinamico.c
	gcc -fopenmp PDSpaceStructure.c libPP_5.0.c EntSalArb_MP.c GNA.c ControlDinamico.c -lfftw3 -lm -o PDStructure.out

PhaseDiagramMP: PhaseDiagramMP.c libPP_5.0.c EntSalArb_MP.c GNA.c ControlDinamico.c
	gcc -fopenmp PhaseDiagramMP.c libPP_5.0.c EntSalArb_MP.c GNA.c ControlDinamico.c -lfftw3 -lm -o PD.out

BrownRemasterOMP: BrownRemasterOMP.c libPP_5.0.c EntSalArb_MP.c GNA.c ControlDinamico.c
	gcc -fopenmp BrownRemasterOMP.c libPP_5.0.c EntSalArb_MP.c GNA.c ControlDinamico.c -lfftw3 -lm -o BrownRem.out
	
CargaYCalculos: CargaYCalculos.c libPP_5.0.c EntSalArb_MP.c GNA.c
	gcc -fopenmp CargaYCalculos.c libPP_5.0.c EntSalArb_MP.c GNA.c -lfftw3 -lm 
	
Competencia: Competencia.c libPP_6.1.c EntSalArb_MP_Comp.c GNA.c
	gcc -fopenmp Competencia.c libPP_6.1.c EntSalArb_MP_Comp.c GNA.c -lfftw3 -lm -o crece.out

BrownRemaster: BrownRemaster.c libPP_4.0.c EntSalArb.c GNA.c
	gcc BrownRemaster.c libPP_4.0.c EntSalArb.c GNA.c -lm

Brown: Brown.c libPP_3.0.c EntSalArb.c GNA.c
	gcc Brown.c libPP_3.0.c EntSalArb.c GNA.c -lm 

CargaTabla: CargaTabla.c libPP_5.0.c EntSalArb_MP.c GNA.c
	gcc CargaTabla.c libPP_5.0.c EntSalArb_MP.c GNA.c -lfftw3 -lm -o CargaTabla.out
