#define K(individual, Model) Model->coagulation_factor * (float)(individual.size)
	
#define R(individual, Model) Model->coagulation_radio_factor*sqrtf(sqrtf((float)(individual.size)))

#define M(individual, Model) Model->metabolic_factor*(float)(individual.size)
