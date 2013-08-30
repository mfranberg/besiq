#ifndef __SNP_COUNT_H__
#define __SNP_COUNT_H__

#include <armadillo>

#include <snp_row.hpp>

/**
 * Counts the number of cases and controls with each genotype. The
 * counts are based on the weight, so an individual with weight 0.5
 * will be counted as 0.5 instead of 1.0.
 *
 * @param row1 The first snp.
 * @param row2 The second snp.
 * @param phenotype The phenotype 0.0 or 1.0.
 * @param weight The weight of each individual.
 * 
 * @return Counts for each genotype. They are represented as a 9x2 matrix,
 *         so that each row is a cell ordered from left to right and top to bottom,
 *         and the first cell in each row is for controls and the second for cases.
 */
arma::mat joint_count(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

/**
 * Counts the number of individuals with each genotype.
 *
 * @param row1 The first snp.
 * @param row2 The second snp.
 * 
 * @return Counts for each genotype. They are represented as a 9 element vector,
 *         so that each cell is ordered from left to right and top to bottom.
 */
arma::vec joint_count(const snp_row &row1, const snp_row &row2);

/**
 * Counts the number of cases and controls with each phenotype. The
 * counts are based on the weight, so an individual with weight 0.5
 * will be counted as 0.5 instead of 1.0.
 *
 * @param row1 The first snp.
 * @param row2 The second snp.
 * @param phenotype The phenotype 0.0 or 1.0.
 * @param weight The weight of each individual.
 * 
 * @return Counts for each phenotype. The are represented as a 2 element vector,
 *         so that the first element is the number of controls and the second the
 *         number of cases. Individuals with the first or second snp missing is
 *         ignored.
 */
arma::vec pheno_count(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

/**
 * Counts the number of cases and controls with each genotype for one
 * of the snps. The counts are based on the weight, so an individual with weight 0.5
 * will be counted as 0.5 instead of 1.0.
 *
 * @param row1 The first snp.
 * @param row2 The second snp.
 * @param phenotype The phenotype 0.0 or 1.0.
 * @param weight The weight of each individual.
 * 
 * @return Counts for each genotype for one snp. They are represented as a 3x2 matrix,
 *         so that each row is a allele (0, 1, 2) and the first cell in each row is 
 *         for controls and the second for cases. Individuals with any snp missing is
 *         ignored.
 */
arma::mat single_count(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

/**
 * Estimates the minor allele frequency for a snp.
 * 
 * @param row A snp.
 * 
 * @return A 3-element vector that contains the probability
 *         of each allele, i.e. it sums to 1.0.
 */
arma::vec compute_maf(const snp_row &row);

#endif /* End of __SNP_COUNT_H__ */
