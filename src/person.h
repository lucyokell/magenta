//
//  MAGENTA
//  person.h
//
//  Created: OJ on 12/01/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing an individidual.
//
// ---------------------------------------------------------------------------


#ifndef PERSON_H
#define PERSON_H

#include <vector>
#include <queue>
#include <deque>
#include "strain.h"
#include "parameters.h"
#include "probability.h"
using namespace std;

class Person {

public:
	enum InfectionStatus
	{
		SUSCEPTIBLE,	// 0
		DISEASED,		// 1
		ASYMPTOMATIC,	// 2
		SUBPATENT,		// 3
		TREATED,		// 4
		PROPHYLAXIS,	// 5
		NUMBER_OF_STATES = 6
	};

	// Person's individual relative biting rate (zeta) due to mosquito biting heterogeneity following a lognormal distribution
	const double m_individual_biting_rate;  // Set for life and thus constant and fine to be public
	const int m_person_ID;					// member variable ID	- fine to be public as constant
	
private:

	static int s_person_ID_generator;		// static member variable for initialising unique identifiers
	int m_person_age;						// Person's age
	InfectionStatus m_infection_state;		// Infection Status enum

	// Person's age dependent biting rate (psi) - See Griffin 2010 S1 for this specific origin
	double m_age_dependent_biting_rate;
	// Person's symptom success rate (phi) - See Griffin 2010 S1 for this specific origin
	double m_symptom_success_rate;

	// Person's transition probabilities, i.e. the probabilty, given their acquired immunity and the treatment seeking rate, that they 
	// become diseased, treated or asymptomatic
	std::vector<double> m_transition_probabilities = {0.33, 0.33, 0.33};

	// Infection transition options. TODO: This does not need to be in every instance of a person and could be declared once globally
	std::vector<InfectionStatus> m_transition_vector { DISEASED, TREATED, ASYMPTOMATIC };

	// Sum of transition probabilities. Will always be 1 unless transitioning from DISEASED, as seems incorrect to transition to ASYMPTOMATIC due to an infection.
	double m_sum_transition_probabilities = 1;

	// Number of bites individual is receiving in this time step - only calculated for those being biten
	int m_number_of_bites = 0;

	// Number of succesful bites individual is receiving in this time step - only calculated for those being biten
	int m_number_of_succesful_bites = 0;

	// Number of strains an individual posseses - this can include multiple same strains
	int m_number_of_strains = 0;

	// Mosquito ID pointers so that we know what mosquitos bit each human
	// TODO: Add this in
	// Mosquito* m_mosquito_pointers;

	// Vector of strain pointers
	std::vector<Strain*> m_active_strain_pointers;

	// Vector of associated active strain acquisition dates, i.e. so we know possible concentrations  of a strain when human is bitten
	std::vector<int> m_active_strain_acquisition_day;

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// IMMUNITY VARIABLES 
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Person's pre-erythrocytic immunity (IB), i.e. infection blocking
	double m_IB;								// infection-blocking immunity 
	double m_biting_success_rate = 0;			// Temporary biting success rate (b) dependent on age and pre-erythocytic immunity
	int m_IB_last_boost_time = 0;				// time that IB was last boosted
	int m_IB_last_calculated_time = 0;			// time IB was last calculated

	// Person's acquired disease blocking immunity 
	double m_ICM;								// maternally acquired immunty (ICM)
	double m_ICA;								// exposure acquired immunity (ICA)
	int m_ICA_last_boost_time = 0;				// time that ICA was last boosted
	int m_I_C_D_CM_last_calculated_time = 0;		// time IC and ID was last calculated

	// Person's blood stage immunity (I.D), i.e. reducing probability of detectiona ond onwards contribution to infectiousness
	double m_ID;								// infection-blocking immunity 
	int m_ID_last_boost_time = 0;				// time that ID was last boosted

	// Person's float for immunities and age to give non integer comparisons for boosting, i.e. so parameter
	// for boosting equals 7.2, it doesn't essentially become 8 for everyone
	double m_immunity_boost_float = runif1(0.0, 1.0);

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// EVENT VARIABLES 
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	int m_day_of_InfectionStatus_change = 0;	// Person's time of InfectionStatus state change
	int m_day_of_strain_clearance = 0;			// Person's time of clearing a strain
	int m_day_of_death = 0;						// Person's time of death
	int m_day_of_next_event = 0;				// Person's closest event day

	std::queue<int> m_infection_time_realisation_queue{};				// First in first out queue to handle multiple infections times that have not been realised yet
	std::queue<InfectionStatus> m_infection_state_realisation_queue{};	// First in first out queue to handle multiple infections states that have not been realised yet
	std::queue<Strain*> m_infection_strain_realisation_queue{};	// First in first out queue to handle multiple infections states that have not been realised yet

	// Variable that allows a check for empty queues when dealing with more than one infection realisation on a day
	int m_infection_realisation_empty_catch = 1;

	

public:

	// Default class constructor
	Person();

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// GETTERS
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	// Get person's age
	int get_m_person_age() { return(m_person_age); }

	// Get person's infection status
	InfectionStatus get_m_infection_state() { return(m_infection_state); }
	
	// Get person's age-dependent immunity
	double get_m_age_dependent_biting_rate() { return(m_age_dependent_biting_rate); }

	// Get person's IB
	double get_m_IB() { return (m_IB); }
	
	// Get person's ICA
	double get_m_ICA() { return (m_ICA); }
	
	// Get person's ICM
	double get_m_ICM() { return (m_ICM); }
	
	// Get person's ID
	double get_m_ID() { return (m_ID); }
	
	// Get person's IB_last_boost_time
	int get_m_IB_last_boost_time() { return (m_IB_last_boost_time); }
	
	// Get person's ICA_last_boost_time
	int get_m_ICA_last_boost_time() { return (m_ICA_last_boost_time); }
	
	// Get person's ID_last_boost_time
	int get_m_ID_last_boost_time() { return (m_ID_last_boost_time); }
	
	// Get person's IB_last_calculated_time
	int get_m_IB_last_calculated_time() { return (m_IB_last_calculated_time); }
	
	// Get person's I_C_D_CM_last_calculated_time
	int get_m_I_C_D_CM_last_calculated_time() { return (m_I_C_D_CM_last_calculated_time); }
	
	// Get person's immunity boost float
	double get_m_immunity_boost_float() { return (m_immunity_boost_float); }
	
	// Get person's time of InfectionStatus state change
	int get_m_day_of_InfectionStatus_change() { return (m_day_of_InfectionStatus_change); }	
	
	// Get person's time of strain clearance
	int get_m_day_of_strain_clearance() { return (m_day_of_strain_clearance); }
	
	// Get person's time of death
	int get_m_day_of_death() { return (m_day_of_death); }
	
	// Get person's infection time realisation queue
	std::queue<int> get_m_infection_time_realisation_queue() { return (m_infection_time_realisation_queue); }	
	
	// Get person's infection state realisation queue
	std::queue<InfectionStatus> get_m_infection_state_realisation_queue() {return (m_infection_state_realisation_queue); }

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// SETTERS
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Set static person generator
	void set_s_person_ID_generator(int x = 1) { Person::s_person_ID_generator = x; }

	// Set person's infection state
	void set_m_infection_state(InfectionStatus x) { m_infection_state = x; }

	// TODO: THESE ARE NOT FORMAL SETTERS AND SHOULD BE REWORDED
	// --------------------------------------------------------------------------------------------------------
	// Set person's individual biting rate
	double set_m_individual_biting_rate(double zeta_meanlog, double zeta_sdlog);

	// Set person's relative chance of biting (psi)
	double set_m_age_dependent_biting_rate(double rho, double a0);

	// Set person's intial age given the maximum and average age
	int set_initial_m_person_age(double average_age);

	// Set person's initial death day
	void set_initial_m_day_of_death(const Parameters &parameters);

	// Set person's transition probabilities
	void set_m_transition_probabilities(double pD, double pT, double pA) { m_transition_probabilities = { pD, pT, pA }; }

	// Set day of next event
	void set_m_day_of_next_event(); 

	// --------------------------------------------------------------------------------------------------------

	// Set person's number of bites being received
	void set_m_number_of_bites(int x) { m_number_of_bites = x; }

	// Set person's number of bites being received
	void set_m_number_of_strains(int x) { m_number_of_strains = x; }

	// Set person's IB
	void set_m_IB(double x) { m_IB = x; }

	// Set person's ID
	void set_m_ID(double x) { m_ID = x; }

	// Set person's ICA
	void set_m_ICA(double x) { m_ICA = x; }

	// Set person's ICM
	void set_m_ICM(double x) { m_ICM = x; }

	// Set person's ICM
	void set_m_day_of_strain_clearance(int x) { m_day_of_strain_clearance = x; }

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// ALLOCATIONS
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Allocate bite to person
	void allocate_bite(const Parameters &parameters);

	// Allocate force of infection to person
	// void allocate_force_of_infection(const Parameters &parameters);

	// Allocate an infection to person, i.e. individual's who return >0 from allocate_force_of_infection()
	void allocate_infection(const Parameters &parameters);

	// Allocate a strain to individual
	void allocate_strain(const Parameters &parameters, Strain* strain_ptr);

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// SCHEDULERS 
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Schedule person's infection state change
	void schedule_m_day_of_InfectionStatus_change(const Parameters &parameters);

	// Schedule person's death day
	void schedule_m_day_of_death(const Parameters &parameters);

	// Schedule person's next strain clearance
	void schedule_m_day_of_strain_clearance(const Parameters &parameters);

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// UPDATERS - functions that happen at the end of a day, i.e. ageing, dying, clearing strains etc. 
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Clear one strain from an individual
	void individual_strain_clearance();

	// Clear all strains from an individual
	void all_strain_clearance();

	// Recover to being susceptible, i.e. clearing all infections and strains and associated timings
	void recover(const Parameters &parameter);

	// Kill person, i.e. reset age to 0, infections to 0, state to susceptible, immunities reset etc
	void die(const Parameters &parameters);

	// Daily update function, i.e. increase ages, set maternal booleans, calculate biting proportions etc
	double update(Parameters &parameters);

	// Event handle, i.e. if any of person's death, state change, strain clearance or state change realisation days are today, respond accordingly
	void event_handle(const Parameters &parameters);
	
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// LOGGERS - functions that report summaries on an individual
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	// Caluclate daily incidence for population and 0-5 years, i.e. would they cause an incident case today.
	// Returns 2 if 0-5 and incident, 1 if just incident, and 0 if neither
	int log_daily_incidence(const Parameters &parameters);

};

#endif // PERSON_H