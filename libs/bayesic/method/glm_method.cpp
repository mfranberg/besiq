#include <bayesic/method/glm_method.hpp>

glm_method::glm_method(method_data_ptr data, const glm_model &model)
: method_type::method_type( data ),
  m_design_matrix( data->phenotype.n_elem, data->covariate_matrix.n_cols + 4 ),
  m_model( model )
{
    /*
     * First three columns are snp1, snp2, snp1 x snp2 and intercept.
     */
    m_design_matrix.col( 3 ) = arma::ones<arma::vec>( data->phenotype.n_elem );
    for(int i = 0; i < data->covariate_matrix.n_cols; i++)
    {
        m_design_matrix.col( i + 4 ) = data->covariate_matrix.col( i );
    }
}

void
glm_method::init(std::ostream &output)
{
    output << "Beta" << "\t" << "SE" << "\t" << "P";
}

void glm_method::run(const snp_row &row1, const snp_row &row2, std::ostream &output)
{
    arma::uvec missing = get_data( )->missing;

    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 && missing[ i ] == 0 )
        {
            m_design_matrix( i, 0 ) = row1[ i ];
            m_design_matrix( i, 1 ) = row2[ i ];
            m_design_matrix( i, 2 ) = row1[ i ] * row2[ i ];
        }
        else
        {
            m_design_matrix( i, 0 ) = 0.0;
            m_design_matrix( i, 1 ) = 0.0;
            m_design_matrix( i, 2 ) = 0.0;
            missing[ i ] = 1;
        }
    }

    irls_info info;
    arma::vec b = irls( m_design_matrix, get_data( )->phenotype, missing, m_model, info );   

    if( info.converged )
    {
        output << b[ 2 ] << "\t" << info.se_beta[ 2 ] << "\t" << info.p_value[ 2 ];
    }
    else
    {
        output << "NA\tNA\tNA";
    }
}
