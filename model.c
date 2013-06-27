
#ifdef SOI
	#define K(individual, Model) (Model->coagulation_factor)
#else
	#define K(individual, Model) (Model->coagulation_factor * (double)((individual.radio)*(individual.radio)))
#endif

#define R(individual, Model) (Model->coagulation_radio_factor*sqrt(sqrt((double)(individual.size))))

#define M(individual, Model) (Model->metabolic_factor*(double)(individual.size))


