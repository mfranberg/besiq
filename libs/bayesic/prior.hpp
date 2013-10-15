#ifndef __PRIOR_H__
#define __PRIOR_H__

#include <vector>

#include <armadillo>

#include <plink/snp_row.hpp>

/**
 * Estimates the prior parameters by first sampling:
 *   1. Taking random pairs of snps.
 *   2. Permuting the phenotype.
 *   3. Taking a random cell.
 *   4. Estimating the risk.
 * 
 * Then over these parameter samples it will fit a beta prior
 * distribution, and return its parameters.
 *
 * @param genotype_matrix Matrix of all genotypes.
 * @param phenotypes The phenotype.
 * @param missing Individuals that will not be included in the analysis,
 *                shuld have missing set to 1.
 * @param num_samples The number of samples to take.
 *
 * @return Suitable beta prior parameters, or 1,1 if they could not
 *         be estimated.
 */
arma::vec estimate_prior_parameters(const std::vector<snp_row> &genotype_matrix, const arma::vec &phenotypes, const arma::uvec &missing, int num_samples);

#endif /* End of __PRIOR_H__ */
