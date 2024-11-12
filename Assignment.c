#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TIME_INTERVAL_FOR_FORCED_DECAY 0.1
#define K_eff 0.25
#define NUMBER_OF_PRECURSOR_GROUPS 10
#define NUETRON_SPEED 22

// -------------------- Estimating Precursors -------------------- //

/*
   -> The first challenge is to sample the decay of precursors in such a way that delayed neutrons are being generated constantly.

   -> Precursor Lifespan can vary from 0.01 sec to 100 sec

   -> At the end of it's lifespan, it decays into a Neutron. This delayed neutron then starts a Neutron Chain, which may form new Precursors in the process

   -> To calculate the average lifetime of the Neutron chain, we will need the probability of a neutron forming a new prompt neutron. 
      This is given by => P_f = K_eff * (1 - Beta), where K_eff is the effective multiplication factor and Beta is the fraction of delayed Neutrons.

   -> If the system is Critical, K_eff = 1
   
   -> Also, if the system is Critical, the probability that the length of Neutron Chain is N is given by P(N) = (P_f ^ (N-1)) * (1 - P_f)
   
   -> The average length of the neutron chain can be given by n_bar = Summation (n * P(n)), where n ranges from 1 to infinity. 
      Solving this gives n_bar = 1 / (1 - P_f)
   
   -> This can lead to large variance (Explained in the paper), and to fix this, we force Precursors to Decay every time interval, leading to new Neutron Chains
   
   -> The probability of a Precursor decaying in time interval t and t + del_t is given by uniform distribution, or p(t) = 1.00 / del_t

   -> For a delay group i, the probability of a precursor having natural decay is given by p_i(t) = lambda_i * e ^ (-lambda_i * (t - t_0))
      Then, for multiple groups, the probability is the summation over all groups

   -> As this precursor is now Biased, being a weighted sum over all precursors, the weight of the delayed neutron should equal the ration of forced probability and natural probability. i.e. w_n(t) = p(t) / Summ(p_i(t)) over all i 
*/

double* precursor_sampling(double* beta_array, double* lambda_array, double t, double t_0) {
   double beta = 0.0;
   double P_t = 0.0;
   for(int i=0; i<NUMBER_OF_PRECURSOR_GROUPS; i++) {
      beta += beta_array[i];
      P_t += beta_array[i] * lambda_array[i] * exp(-1.00 * (t - t_0) * lambda_array[i]);
   }
   P_t /= beta;

   double P_f = K_eff * (1.00 - beta);
   double n_bar = 1.00 / (1.00 - P_f);
   double P_bar = 1.00 / TIME_INTERVAL_FOR_FORCED_DECAY;
   double w_n = P_t / P_bar;

   double* arr = malloc(5 * sizeof(double));
   arr[0] = P_f;
   arr[1] = n_bar;
   arr[2] = P_bar;
   arr[3] = P_t;
   arr[4] = w_n;

   return arr;
}


// -------------------- Precursors Spatial Distribution -------------------- //

/*
   -> When initializing a system, first a criticality calculation is done until the source has converged. This is the steady state neutron distribution and from this distribution the precursor and prompt neutron distribution can be calculated.

   -> As we are assumiung all the precursors to be in the same delay group, 
      C_0(r) = (beta / lambda_b) neu * Sigma_f * phi(x), where lambda_b is the averaged lambda.

   -> Therefore, fraction of prompt neutrons at a particular r is given by - 1 / (1 + (beta / lambda) * v * neu * Sigma_f)

*/

double precursor_spatial_distribution(double neu, double sigma_f, double* beta_array, double* lambda_array) {
   double lambda_b = 0.0;
   double beta = 0.0;
   double sum_of_beta_by_lambda = 0.0;
   for(int i=0; i<NUMBER_OF_PRECURSOR_GROUPS; i++) {
      beta += beta_array[i];
      sum_of_beta_by_lambda += beta_array[i] / lambda_array[i];
   }
   lambda_b = beta / sum_of_beta_by_lambda;

   return 1.00 / (1 + (beta / lambda_b) * NUETRON_SPEED * neu * sigma_f);
}


// -------------------- Precursors Time Distribution -------------------- //

/*
   -> With exponential decay a particle has no age. There is always the same probability that a particle has its next decay, no matter what time it has lived before. In this case however as stated before a combined precursor particle does not have pure exponential decay.

   -> The importance of a precursor group varies over time, and is given by 
      C_i(t) / C(t) = ((beta_i / beta) * e ^ (-lambda_i * t)) / [Summation of ((beta_j / beta) * e ^ (-lambda_j * t)) over all j]

   -> At a stationary condition, the fraction of precursors per group remains the same, and can be given by 
      C_i0(t) / C_0(t) = (beta_i * lambda ^ b) / (beta * lambda_i)

   -> To achieve this in a Monte Carlo simulation the combined precursor particles should be started with an age between −∞ and 0. This way the ratios between the different precursor groups are correct.
*/

double* precursor_time_distribution_at_time_t(double* lambda_array, double* beta_array, double t, double t_0) {
   double* importance_of_precursor_group = malloc(sizeof(double) * NUMBER_OF_PRECURSOR_GROUPS);
   double beta = 0.0;
   double C_at_t = 0.0; 
   for(int i = 0; i<NUMBER_OF_PRECURSOR_GROUPS; i++) {
      beta += beta_array[i];
      C_at_t += beta_array[i] * exp(-1.00 * lambda_array[i] * t);
   }
   
   C_at_t /= beta;
   
   for(int i = 0; i<NUMBER_OF_PRECURSOR_GROUPS; i++) {
      importance_of_precursor_group[i] = ((beta_array[i] / beta) * exp(-1.00 * lambda_array[i] * t)) / C_at_t;
   }
   return importance_of_precursor_group;
}

double* precursors_per_group_at_stationary_case(double* lambda_array, double* beta_array) {
   double* stationary_case_fraction_of_precursor_group = malloc(sizeof(double) * NUMBER_OF_PRECURSOR_GROUPS);
   double beta = 0.0;
   double lambda_b = 0.0;
   double sum_of_beta_by_lambda = 0;
   for(int i = 0; i<NUMBER_OF_PRECURSOR_GROUPS; i++) {
      beta += beta_array[i];
      sum_of_beta_by_lambda += beta_array[i] / lambda_array[i];
   }
   lambda_b = beta / sum_of_beta_by_lambda;
   
   for(int i=0; i<NUMBER_OF_PRECURSOR_GROUPS; i++) {
      stationary_case_fraction_of_precursor_group[i] = (beta_array[i] * lambda_b) / (beta * lambda_array[i]);
   }
   return stationary_case_fraction_of_precursor_group;
}

double calculation_of_effective_weight(double t_0, double* lambda_array, double* beta_array) {
   double sum = 0.0;
   double beta = 0;

   for(int i=0; i<NUMBER_OF_PRECURSOR_GROUPS; i++) {
      sum += beta_array[i] * exp(lambda_array[i] * t_0);
      beta += beta_array[i];
   }
   return sum / beta;
}


// -------------------- Calculation Scheme -------------------- //

/*
   -> To ensure that the Monte Carlo method can adapt to its changing environment, the simulation is split in time intervals. The size of these time steps can be chosen freely, but all precursors are forced to have a decay in every time step. If a particle crosses the time boundary of a time step, it is stored for the next time interval.

   -> After each time interval the system properties can be adjusted

   ->
*/

void calculate(double simulation_time, double t_0, double* beta_array, double* lambda_array, int* partitions, long long number_of_particles) {
   double beta = 0;
   for(int i=0; i<NUMBER_OF_PRECURSOR_GROUPS; i++) {
      beta += beta_array[i];
   }
   int initial_neutron_count = 0;
   double prob_neutron_produces_prompt_neutron = (1 - beta) * K_eff;

   for(int t=t_0; t<simulation_time; t += TIME_INTERVAL_FOR_FORCED_DECAY) {
      for(int i=0; i<number_of_particles; i++) {
         for(int j=0; j<NUMBER_OF_PRECURSOR_GROUPS; j++) {
            if(partitions[j] < i && i <= partitions[j+1]) {
               double decay_probability = lambda_array[j] * exp(-1.00 * lambda_array[j] * (t - t_0));
               double random_number_for_decay = ((double)rand()) / RAND_MAX;
               double random_number_for_prompt_neutron_generation = ((double)rand()) / RAND_MAX;
               if(random_number_for_decay < decay_probability) initial_neutron_count++;
               if(random_number_for_prompt_neutron_generation < prob_neutron_produces_prompt_neutron) {
                  partitions[NUMBER_OF_PRECURSOR_GROUPS+1]++;
                  number_of_particles++;   
               }
            } 
         }
      }
   }
}


int main(void) {
   srand(time(NULL));
}