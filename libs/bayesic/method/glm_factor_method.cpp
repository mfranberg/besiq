#include <bayesic/method/glm_factor_method.hpp>

#include <dcdflib/libdcdf.hpp>

glm_factor_method::glm_factor_method(method_data_ptr data, const glm_model &model)
: method_type::method_type( data ),
  m_null_design_matrix( data->phenotype.n_elem, data->covariate_matrix.n_cols + 5 ),
  m_alt_design_matrix( data->phenotype.n_elem, data->covariate_matrix.n_cols + 9 ),
  m_model( model )
{
    /*
     * Null matrix.
     */
    m_null_design_matrix.col( 4 ) = arma::ones<arma::vec>( data->phenotype.n_elem );
    for(int i = 0; i < data->covariate_matrix.n_cols; i++)
    {
        m_null_design_matrix.col( i + 5 ) = data->covariate_matrix.col( i );
    }

    /*
     * Alternative matrix.
     */
    m_alt_design_matrix.col( 8 ) = arma::ones<arma::vec>( data->phenotype.n_elem );
    for(int i = 0; i < data->covariate_matrix.n_cols; i++)
    {
        m_alt_design_matrix.col( i + 9 ) = data->covariate_matrix.col( i );
    }
}

void
glm_factor_method::init(std::ostream &output)
{
    output << "LR" << "\t" << "P";
}

void
glm_factor_method::init_alt(const snp_row &row1, const snp_row &row2, arma::mat &design_matrix, arma::uvec &missing)
{
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 && missing[ i ] == 0 )
        {
            double s11 = row1[ i ] == 1 ? 1.0 : 0.0;
            double s12 = row1[ i ] == 2 ? 1.0 : 0.0;
            double s21 = row2[ i ] == 1 ? 1.0 : 0.0;
            double s22 = row2[ i ] == 2 ? 1.0 : 0.0;

            design_matrix( i, 0 ) = s11;
            design_matrix( i, 1 ) = s12;
            design_matrix( i, 2 ) = s21;
            design_matrix( i, 3 ) = s22;
            design_matrix( i, 4 ) = s11*s21;
            design_matrix( i, 5 ) = s11*s22;
            design_matrix( i, 6 ) = s12*s21;
            design_matrix( i, 7 ) = s12*s22;
        }
        else
        {
            design_matrix( i, 0 ) = 0.0;
            design_matrix( i, 1 ) = 0.0;
            design_matrix( i, 2 ) = 0.0;
            design_matrix( i, 3 ) = 0.0;
            design_matrix( i, 4 ) = 0.0;
            design_matrix( i, 5 ) = 0.0;
            design_matrix( i, 6 ) = 0.0;
            design_matrix( i, 7 ) = 0.0;
            missing[ i ] = 1;
        }
    }
}

void
glm_factor_method::init_null(const snp_row &row1, const snp_row &row2, arma::mat &design_matrix, arma::uvec &missing)
{
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 && missing[ i ] == 0 )
        {
            double s11 = row1[ i ] == 1 ? 1.0 : 0.0;
            double s12 = row1[ i ] == 2 ? 1.0 : 0.0;
            double s21 = row2[ i ] == 1 ? 1.0 : 0.0;
            double s22 = row2[ i ] == 2 ? 1.0 : 0.0;

            design_matrix( i, 0 ) = s11;
            design_matrix( i, 1 ) = s12;
            design_matrix( i, 2 ) = s21;
            design_matrix( i, 3 ) = s22;
        }
        else
        {
            design_matrix( i, 0 ) = 0.0;
            design_matrix( i, 1 ) = 0.0;
            design_matrix( i, 2 ) = 0.0;
            design_matrix( i, 3 ) = 0.0;
            missing[ i ] = 1;
        }
    }
}

void glm_factor_method::run(const snp_row &row1, const snp_row &row2, std::ostream &output)
{
    arma::uvec missing = get_data( )->missing;

    irls_info null_info;
    init_null( row1, row2, m_null_design_matrix, missing );
    irls( m_null_design_matrix, get_data( )->phenotype, missing, m_model, null_info );

    irls_info alt_info;
    init_alt( row1, row2, m_alt_design_matrix, missing );
    irls( m_alt_design_matrix, get_data( )->phenotype, missing, m_model, alt_info );

    if( null_info.converged && alt_info.converged )
    {
        double LR = -2 * ( null_info.logl - alt_info.logl );
        double p = 1.0 - chi_square_cdf( LR, 4 );

        output << LR << "\t" << p;
    }
    else
    {
        output << "NA\tNA";
    }
}
