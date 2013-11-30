
#ifdef SOI
	#define K(individual, Model) (Model->coagulation_factor)
#else
	#define K(individual, Model) (Model->coagulation_factor * (double)((individual.radio)*(individual.radio)))
#endif

//#define R(individual, Model) (Model->coagulation_radio_factor*pow((double)(individual.size), 0.375))
//#define M(individual, Model) (Model->metabolic_factor*(double)(individual.size))

#define R(individual, Model) (individual.size/2.0)
#define M(individual, Model) (Model->metabolic_factor*pow((double)(individual.size),2.6666))



