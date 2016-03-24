#ifndef __COVARIATE_H__
#define __COVARIATE_H__

#include <iostream>

#include <armadillo>

struct pio_sample_t;

/**
 * Parses the contents of a csv stream and returns them
 * as a matrix. Missing values will be set to 1 in the
 * given vector.
 *
 * @param stream The stream for the csv data.
 * @param missing The individuals with missing values will be set to 1,
 *                others will remain untouched.
 * @param order This vector defines the order of the individuals that will
 *              be parsed from the covariate file.
 * @param out_header The header names will be stored here.
 * @param missing_string The string that indicates a missing value.
 *
 * @return A matrix that contains the parsed covariates.
 */
arma::mat
parse_covariate_matrix(std::istream &stream, arma::uvec &missing, const std::vector<std::string> &order, std::vector<std::string> *out_header = NULL, const char *missing_string = "NA");

/**
 * Parsers phenotypes from a csv stream and returns them as a 
 * a matrix. Missing values will be set to 1 in the given vector.
 *
 * @param stream The stream to read phenotypes from.
 * @param missing The individuals with missing values will be set to 1,
 *                others will remain untouched.
 * @param order This vector defines the order of the individuals that will
 *              be parsed from the phenotype file.
 * @param pheno_name If different than "" it will search for the column with this
 *                   name and return it.
 * @param missing_string The string that indicates a missing value.
 *
 * @return A vector containing the parsed phenotypes.
 */
arma::vec
parse_phenotypes(std::istream &stream, arma::uvec &missing, const std::vector<std::string> &order, std::string pheno_name = "", const char *missing_string = "NA");

/**
 * Parsers environment variables from a csv stream and returns them as a 
 * a matrix. Missing values will be set to 1 in the given vector.
 *
 * @param stream The stream to read phenotypes from.
 * @param missing The individuals with missing values will be set to 1,
 *                others will remain untouched.
 * @param order This vector defines the order of the individuals that will
 *              be parsed from the environment file.
 * @param env_name If different than "" it will search for the column with this
 *                   name and return it.
 * @param missing_string The string that indicates a missing value.
 *
 * @return A matrix containing the parsed environmental variables.
 */
arma::mat
parse_env(std::istream &stream, arma::uvec &missing, const std::vector<std::string> &order, std::vector<std::string> *out_header = NULL, std::string env_name = "", const char *missing_string = "NA");

/**
 * Parses an environmental factor from a csv stream and returns them as a 
 * a matrix. Missing values will be set to 1 in the given vector.
 *
 * @param stream The stream to read the environmental factor from.
 * @param missing The individuals with missing values will be set to 1,
 *                others will remain untouched.
 * @param order This vector defines the order of the individuals that will
 *              be parsed from the phenotype file.
 * @param levels The number of levels of the factor, 1 is continuous, 2
 *               is binary etc.
 * @param missing_string The string that indicates a missing value.
 *
 * @return A vector containing the parsed environmental variables.
 */
arma::mat
parse_environment(std::istream &stream, arma::uvec &missing, const std::vector<std::string> &order, unsigned int levels, const char *missing_string = "NA");

/**
 * Creates a vector of phenotypes from the given sample vector.
 *
 * @param samples Plink samples.
 * @param missing Missing values will be added here.
 *
 * @return A vector of phenotypes of the same length as the input samples.
 */
arma::vec create_phenotype_vector(const std::vector<pio_sample_t> &samples, arma::uvec &missing);

#endif /* End of __COVARIATE_H__ */
