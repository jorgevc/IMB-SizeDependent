#define K(individual, Model) Model->coagulation_factor * ((float)(individual.size) + (Model->health_factor)*(float)(individual.size*individual.health))
	
#define R(individual, Model) Model->coagulation_radio_factor*sqrtf((float)(individual.size))

#define M(individual, Model) Model->metabolic_factor*(pow((float)(individual.size + individual.health), 2.0))
