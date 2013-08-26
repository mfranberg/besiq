#ifndef __SNP_COUNT_H__
#define __SNP_COUNT_H__

#include <armadillo>

#include <plink_file.hpp>

/**
 * Counts the number cases and controls with each genotype.
 */
arma::mat joint_count(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

arma::vec joint_count(const snp_row &row1, const snp_row &row2);

arma::vec pheno_count(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

arma::mat single_count(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

arma::vec compute_maf(const snp_row &row);

#endif /* End of __SNP_COUNT_H__ */
