#include <vector>

#include <armadillo>

#include <plink/plink_file.hpp>

/**
 * This class is responsible for both genetic and environmental data, along
 * with computing certain properties of this data that relates to the LARS
 * algorithm. Primarily, to avoid passing and operating on each data type
 * separately within the algorithm.
 */
class gene_environment
{
public:
    /**
     * Constructor.
     *
     * @param genotypes A matrix of genotypes.
     * @param cov A matrix of covariates.
     * @param phenotype A vector of phenotypes.
     * @param cov_names Names of the covariates.
     */
    gene_environment(genotype_matrix_ptr genotypes, const arma::mat &cov, const arma::vec &phenotype, const std::vector<std::string> &cov_names);
    /**
     * Imputes the missing genotypes, covariates and phenotypes.
     */
    void impute_missing();

    /**
     * Centers a phenotype and returns it.
     *
     * @return Returns a centered phenotype.
     */
    arma::vec get_centered_phenotype();

    /**
     * Returns the number of samples.
     *
     * @return the number of samples.
     */
    size_t get_num_samples();

    /**
     * Returns the number of variables.
     *
     * @return the number of variables.
     */
    size_t get_num_variables();
   
    /**
     * Returns the name of the variable with the given index.
     *
     * @param index Index of the variable whoose name to return.
     *
     * @return The name of the variable with the given index.
     */ 
    std::string get_name(size_t index);
   
    /**
     * Returns the complete list of variable names.
     *
     * @return List of variables.
     */ 
    std::vector<std::string> get_names();

    /**
     * Calculates the correlation between each variable and
     * the given vector of residuals. The results are stored
     * in the supplied vector c.
     *
     * @param residual A vector of residuals.
     * @param c A vector of at least the size of the number of variables,
     *          to store the computed correlations in.
     */
    void calculate_cor(const arma::vec &residual, arma::vec &c);
    
    /**
     * Computes the vector a from the LARS paper from the
     * vector u.
     *
     * @param u The u-vector from the LARS paper.
     *
     * @return The a-vector from the LARS paper.
     */
    arma::vec eig_prod(const arma::vec &u);

    /**
     * Returns a matrix of standardized variables according to the
     * given vector of indicides of active variables.
     *
     * @param active  Indicies of active variables.
     *
     * @return Matrix of standardized active variables.
     */
    arma::mat get_active(const arma::uvec &active);
    
    /**
     * Returns a matrix of raw variables according to the
     * given vector of indicides of active variables.
     *
     * @param active  Indicies of active variables.
     *
     * @return Matrix of standardized active variables.
     */
    arma::mat get_active_raw(const arma::uvec &active);

private:
    /**
     * Assigns missing genotypes with genotype mean.
     */
    void fill_missing_genotypes();

    /**
     * Assigns missing phenotypes with phenotype mean.
     */
    void fill_missing_phenotypes();

    /**
     * Assigns missing covariates with covariate mean.
     */
    void fill_missing_cov();

    /**
     * Computes the mean and standard deviation of each variable
     * and stores it in the m_mean and m_sd vectors. This is performed
     * in this way because it is too expensive to store the genotypes
     * as floats.
     */
    void compute_mean_sd();

private:
    /**
     * Matrix of genotypes.
     */
    genotype_matrix_ptr m_genotypes;

    /**
     * Matrix of covariates.
     */
    arma::mat m_cov;

    /**
     * Vector of phenotypes.
     */
    arma::vec m_phenotype;

    /**
     * Vector of mean values for each variable.
     */
    arma::vec m_mean;

    /**
     * Vector of standard deviation for each variable.
     */
    arma::vec m_sd;

    /**
     * Vector variable names.
     */
    std::vector<std::string> m_names;
};
