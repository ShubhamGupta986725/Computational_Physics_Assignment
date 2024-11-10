#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TIME_INTERVAL_FOR_FORCED_DECAY 0.1

// -------------------- NOTES -------------------- //

/*
   -> Assuming all Precursors to be in the same singular Delay Group (i = 1)



*/




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

double* precursor_sampling(double k_eff, double beta, double lambda, double t, double t_0) {
   double P_f = k_eff * (1.00 - beta);
   double n_bar = 1.00 / (1.00 - P_f);
   double P_bar = 1.00 / TIME_INTERVAL_FOR_FORCED_DECAY;
   double P_i = lambda * exp(-1.00 * lambda * (t - t_0));
   double w_n = P_i / P_bar;

   double arr[] = [P_f, n_bar, P_bar, P_i, w_n];
   return arr;
}


// -------------------- Precursors Spatial Distribution -------------------- //

/*

   -> When initializing a system, first a criticality calculation is done until the source has converged. This is the steady state neutron distribution and from this distribution the precursor and prompt neutron distribution can be calculated.

   -> As we are assumiung all the precursors to be in the same delay group, 
      C_0(r) = (beta / lambda_b) neu * Sigma_f * phi(x), where lambda_b is the averaged lambda. As we have a singular delay group, lambda_b = lambda

   -> Therefore, fraction of prompt neutrons at a particular r is given by - 1 / (1 + (beta / lambda) * v * neu * Sigma_f)

*/

void precursor_spatial_distribution(double neu, double v, double sigma_f, double beta, double lambda) {
   return 1.00 / (1 + (beta / lambda) * v * neu * sigma_f);
}


int main(void) {
   printf("%lf\n", TIME_INTERVAL_FOR_FORCED_DECAY);

}