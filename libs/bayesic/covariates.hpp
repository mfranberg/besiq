#ifndef __COVARIATE_H__
#define __COVARIATE_H__

#include <iostream>

#include <armadillo>

/**
 * This function assumes that 
 *
 */
arma::mat
parse_covariate_matrix(std::istream &stream, const char *missing_string = "NA");

#endif /* End of __COVARIATE_H__ */
