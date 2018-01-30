#include <vector>

#include <armadillo>

#include <plink/plink_file.hpp>

class lars_variables
{
public:
    virtual arma::vec get_centered_phenotype( ) const = 0;
    virtual size_t get_num_samples( ) const = 0;
    virtual size_t get_num_variables( ) const = 0;
    virtual std::string get_name(size_t index) const = 0;
    virtual void calculate_cor(const arma::vec &residual, arma::vec &c) const = 0;
    virtual arma::vec eig_prod(const arma::vec &u) const = 0;
    virtual arma::mat get_active(const arma::uvec &active) const = 0;
};

class null_lars : public lars_variables
{
    public:
        null_lars(const arma::mat &X, const arma::vec &y)
            : m_X( X ),
              m_y( y )
        {
        }
    
        arma::vec get_centered_phenotype() const
        {
            return m_y - arma::mean( m_y );
        }
        
        size_t get_num_samples() const
        {
            return m_y.n_elem;
        }

        size_t get_num_variables() const
        {
            return m_X.n_cols;
        }
        
        std::string get_name(size_t index) const
        {
            return "";
        }
    
        void calculate_cor(const arma::vec &residual, arma::vec &c) const
        {
            for(int i = 0; i < m_X.n_cols; i++)
            {
                c[ i ] = arma::dot( m_X.col( i ), residual );
            }
        }

        arma::vec eig_prod(const arma::vec &u) const
        {
            return m_X.t( ) * u;
        }
        
        arma::mat get_active(const arma::uvec &active) const
        {
            return m_X.cols( active );
        }

    private:
        arma::mat m_X;
        arma::vec m_y;
};

/**
 * This class is responsible for both genetic and environmental data, along
 * with computing certain properties of this data that relates to the LARS
 * algorithm. Primarily, to avoid passing and operating on each data type
 * separately within the algorithm.
 */
class gene_environment : public lars_variables
{
public:
    /**
     * Constructor.
     *
     * @param genotypes A matrix of genotypes.
     * @param cov A matrix of covariates.
     * @param phenotype A vector of phenotypes.
     * @param cov_names Names of the covariates.
     * @param only_main Exclude gene-environment.
     */
    gene_environment(genotype_matrix_ptr genotypes, const arma::mat &cov, const arma::vec &phenotype, const std::vector<std::string> &cov_names, bool only_main);
    /**
     * Imputes the missing genotypes, covariates and phenotypes.
     */
    void impute_missing();

    /**
     * Centers a phenotype and returns it.
     *
     * @return Returns a centered phenotype.
     */
    arma::vec get_centered_phenotype() const;

    /**
     * Returns the number of samples.
     *
     * @return the number of samples.
     */
    size_t get_num_samples() const;

    /**
     * Returns the number of variables.
     *
     * @return the number of variables.
     */
    size_t get_num_variables() const;
   
    /**
     * Returns the name of the variable with the given index.
     *
     * @param index Index of the variable whoose name to return.
     *
     * @return The name of the variable with the given index.
     */ 
    std::string get_name(size_t index) const;
   
    /**
     * Returns the complete list of variable names.
     *
     * @return List of variables.
     */ 
    std::vector<std::string> get_names() const;

    /**
     * Returns the number of minor alleles for the given index, and
     * 0 if index is out of bounds.
     *
     * @param index Index to the genotype.
     * 
     * @return The number of minor alleles for the given index, and
     * 0 if the index is out of bounds.
     */
    unsigned int compute_num_minor(size_t index) const;

    /**
     * Calculates the correlation between each variable and
     * the given vector of residuals. The results are stored
     * in the supplied vector c.
     *
     * @param residual A vector of residuals.
     * @param c A vector of at least the size of the number of variables,
     *          to store the computed correlations in.
     */
    void calculate_cor(const arma::vec &residual, arma::vec &c) const;
    
    /**
     * Computes the vector a from the LARS paper from the
     * vector u.
     *
     * @param u The u-vector from the LARS paper.
     *
     * @return The a-vector from the LARS paper.
     */
    arma::vec eig_prod(const arma::vec &u) const;

    /**
     * Returns a matrix of standardized variables according to the
     * given vector of indicides of active variables.
     *
     * @param active  Indicies of active variables.
     *
     * @return Matrix of standardized active variables.
     */
    arma::mat get_active(const arma::uvec &active) const;
    
    /**
     * Returns a matrix of raw variables according to the
     * given vector of indicides of active variables.
     *
     * @param active  Indicies of active variables.
     *
     * @return Matrix of standardized active variables.
     */
    arma::mat get_active_raw(const arma::uvec &active) const;

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

    /**
     * Only consider main effects.
     */
    bool m_only_main;
};
