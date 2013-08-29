#ifndef __COVARIATE_H__
#define __COVARIATE_H__

#include <iostream>

#include <armadillo>

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
 * @param missing_string The string that indicates a missing value.
 *
 * @return A matrix that contains the parsed covariates.
 */
arma::mat
parse_covariate_matrix(std::istream &stream, arma::uvec &missing, const std::vector<std::string> &order, const char *missing_string = "NA");

/**
 * Parsers phenotypes from a csv stream and returns them as a 
 * a matrix. Missing values will be set to 1 in the given vector.
 *
 * @param stream The stream to read phenotypes from.
 * @param missing The individuals with missing values will be set to 1,
 *                others will remain untouched.
 * @param order This vector defines the order of the individuals that will
 *              be parsed from the phenotype file.
 * @param missing_string The string that indicates a missing value.
 *
 * @return A vector containing the parsed phenotypes.
 */
arma::vec
parse_phenotypes(std::istream &stream, arma::uvec &missing, const std::vector<std::string> &order, const char *missing_string = "NA");

#endif /* End of __COVARIATE_H__ */
