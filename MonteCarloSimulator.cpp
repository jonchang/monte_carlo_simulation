#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List SimulateCPBDP(double time,                // the time to run the process for
                         double init_lambda,         // the initial speciation rate
                         double init_mu,             // the initial extinction rate
                         double transition_rate,     // the rate at which new processes arise
                         double lambda_prior,        // the rate of the speciation-rate exponential prior
                         double mu_prior){           // the rate of the extinction-rate exponential prior
  
  double current_time = 0;
  
  // initialize the first process.
  int number_of_processes = 1;
  
  // initialize the vector that tracks the number of species.
  std::vector<int> number_of_species;
  int num_of_species = 1;
  
  // put the first species in the vector.
  number_of_species.push_back(1);
  
  // initialize the speciation and extinction rate vectors, place initial values.
  std::vector<double> speciation_rates;
  speciation_rates.push_back(init_lambda);
  std::vector<double> extinction_rates;
  extinction_rates.push_back(init_mu);

  double return_time = 0;
  
  while (true) {
    
    int which_min_event = 0;
    double next_event_time = INFINITY;
    
    // get the next event time.
    for (int i = 0; i < number_of_processes; ++i) {
      double t = Rcpp::rexp(1,number_of_species[i] * (transition_rate + speciation_rates[i] + extinction_rates[i]))[0];
      if ( t < next_event_time ) {
        which_min_event = i;
        next_event_time = t;
      }
    }
    
    // increment time
    current_time += next_event_time;
    
    // if we have exceeded the time allowed, end the simulation
    if ( current_time >= time ) {
      return_time = time;
      goto ret_label;
    }
    
    // get a random uniform number
    double u = Rcpp::runif(1,0,1)[0];
    
    // get the current total rate for the process that has the event
    double total_rate = transition_rate + speciation_rates[which_min_event] + extinction_rates[which_min_event];
    
    // perform the relevant event on the process
    if ( u < (transition_rate / total_rate) ) { // transition event
      
      number_of_processes++; // increment the number of processes
      number_of_species[which_min_event]--; // remove the species with the new process from its current process
      number_of_species.push_back(1); // and place it in the new process
      double new_speciation_rate = Rcpp::rexp(1,lambda_prior)[0]; // generate the speciation rate for the new process
      double new_extinction_rate = Rcpp::rexp(1,mu_prior)[0]; // generate the extinction rate for the new process
      
      // add the new process
      speciation_rates.push_back(new_speciation_rate); // add the speciation rate for the new process
      extinction_rates.push_back(new_extinction_rate); // add the extinction rate for the new process
      
    } else if ( u < ((transition_rate + speciation_rates[which_min_event]) / total_rate) ) { // speciation event
      
      number_of_species[which_min_event]++; // increment the number of species for this process
      num_of_species++; // increment the total number of species across processes
      
    } else { // extinction event
      
      number_of_species[which_min_event]--; // decement the number of species for this process
      num_of_species--; // decement the total number of species across processes
      
    }
    
    // if the previous process is extinct, remove it
    if ( number_of_species[which_min_event] == 0 ) {
      number_of_processes--; // decrement the number of processes
      number_of_species.erase(number_of_species.begin()+which_min_event); // delete the old process
      speciation_rates.erase(speciation_rates.begin()+which_min_event); // delete the old speciation rate
      extinction_rates.erase(extinction_rates.begin()+which_min_event); // delete the old extinction rate
    }
    
    // if the total number of species is 0, end the simulation and return the current time
    if ( num_of_species == 0 ) {
      return_time = current_time;
      goto ret_label;
    }
    
    // if the total number of species is large, check to see if the probability of 
    // eventual extinction of all the species (under a constant-rate assumption) is
    // very small. if so, terminate the process. this prevents the birth-death
    // process from exploding towards infinity, where the probability of ultimate
    // extinction approaches 0, and computation becomes excessively slow.
    // do this check each time the number of species is divisible by 1000.
    
    if ( num_of_species % 1000 == 0 ) {
      // initialize the probability of ultimate extinction to 1.
      double p = 1;
      for(int i = 0; i < number_of_processes; ++i){ // for each process
        // compute the probability that a _single_ species in that process goes extinct, assuming its rate does not change
        double prob = 1 - (speciation_rates[i] - extinction_rates[i]) / ( speciation_rates[i] - extinction_rates[i] * std::exp( ( extinction_rates[i] - speciation_rates[i] ) * (time - current_time) ) );
        // raise that probability to the number of species in the process, and multiply the 
        // probability of extinction by that number
        p *= pow(prob,number_of_species[i]);
      }
      // if the probability of ultimate extinction of _every_ species is less than 1 in a billion, end the simulation
      if ( p < 1e-9 ) {
        return_time = time;
        goto ret_label;
      }
    }
    
  }

ret_label:
  Rcpp::DataFrame ret_df = Rcpp::DataFrame::create(
          Rcpp::_["number_of_species"]= number_of_species,
          Rcpp::_["speciation_rates"]= speciation_rates,
          Rcpp::_["extinction_rates"]= extinction_rates);
  return Rcpp::List::create(
          Rcpp::_["time"] = return_time,
          Rcpp::_["data"] = ret_df);
}
