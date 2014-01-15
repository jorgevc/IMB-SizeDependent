
#ifdef SOI
	#define K(individual, Model) (Model->coagulation_factor)
#else
	#define K(individual, Model) (Model->coagulation_factor * (double)((individual.radio)*(individual.radio)))
#endif

#define R(individual, Model) (individual.size_float/2.0)
#define M(individual, Model) (Model->metabolic_factor*pow(individual.size_float,2.6666))



