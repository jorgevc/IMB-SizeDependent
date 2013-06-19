
#ifdef SOI
	#define K(individual, Model) (Model->coagulation_factor)
#else
	#define K(individual, Model) (Model->coagulation_factor * (float)((individual.radio)*(individual.radio)))
#endif

#define R(individual, Model) (Model->coagulation_radio_factor*sqrtf(sqrtf((float)(individual.size))))

#define M(individual, Model) (Model->metabolic_factor*(float)(individual.size))


