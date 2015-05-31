#include <besiq/model_matrix.hpp>

general_matrix::general_matrix(const arma::mat &cov, size_t n, size_t num_null, size_t num_alt)
    : m_null( n, cov.n_cols + num_null ),
      m_alt( n, cov.n_cols + num_alt )
{
    /*
     * Null matrix.
     */
    m_null.col( num_null - 1 ) = arma::ones<arma::vec>( n );
    for(int i = 0; i < cov.n_cols; i++)
    {
        m_null.col( i + num_null ) = cov.col( i );
    }

    /*
     * Alternative matrix.
     */
    m_alt.col( num_alt - 1 ) = arma::ones<arma::vec>( n );
    for(int i = 0; i < cov.n_cols; i++)
    {
        m_alt.col( i + num_alt ) = cov.col( i );
    }

    m_num_null = num_null;
    m_num_alt = num_alt;
}

general_matrix::~general_matrix()
{
}

const arma::mat &
general_matrix::get_alt()
{
    return m_alt;
}

const arma::mat &
general_matrix::get_null()
{
    return m_null;
}

size_t
general_matrix::num_df()
{
    return num_alt( ) - num_null( );
}

size_t
general_matrix::num_alt()
{
    return m_num_alt;
}

size_t
general_matrix::num_null()
{
    return m_num_null;
}

additive_matrix::additive_matrix(const arma::mat &cov, size_t n)
    : general_matrix( cov, n, 3, 4 )
{
}

void
additive_matrix::update_matrix(const snp_row &row1, const snp_row &row2, arma::uvec &missing)
{
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 && missing[ i ] == 0 )
        {
            m_alt( i, 0 ) = row1[ i ];
            m_alt( i, 1 ) = row2[ i ];
            m_alt( i, 2 ) = row1[ i ] * row2[ i ];
            
            m_null( i, 0 ) = row1[ i ];
            m_null( i, 1 ) = row2[ i ];
        }
        else
        {
            m_alt( i, 0 ) = 0.0;
            m_alt( i, 1 ) = 0.0;
            m_alt( i, 2 ) = 0.0;
            
            m_null( i, 0 ) = 0.0;
            m_null( i, 1 ) = 0.0;
            
            missing[ i ] = 1;
        }
    }
}

tukey_matrix::tukey_matrix(const arma::mat &cov, size_t n)
    : general_matrix( cov, n, 5, 6 )
{
}

void
tukey_matrix::update_matrix(const snp_row &row1, const snp_row &row2, arma::uvec &missing)
{
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 && missing[ i ] == 0 )
        {
            double s11 = row1[ i ] == 1 ? 1.0 : 0.0;
            double s12 = row1[ i ] == 2 ? 1.0 : 0.0;
            double s21 = row2[ i ] == 1 ? 1.0 : 0.0;
            double s22 = row2[ i ] == 2 ? 1.0 : 0.0;

            m_alt( i, 0 ) = s11;
            m_alt( i, 1 ) = s12;
            m_alt( i, 2 ) = s21;
            m_alt( i, 3 ) = s22;
            // Note: at most one of these terms are 1, so this is either 0 or 1.
            m_alt( i, 4 ) = s11*s21 + s11*s22 + s12*s21 + s12*s22;
            
            m_null( i, 0 ) = s11;
            m_null( i, 1 ) = s12;
            m_null( i, 2 ) = s21;
            m_null( i, 3 ) = s22;
        }
        else
        {
            m_alt( i, 0 ) = 0.0;
            m_alt( i, 1 ) = 0.0;
            m_alt( i, 2 ) = 0.0;
            m_alt( i, 3 ) = 0.0;
            m_alt( i, 4 ) = 0.0;
            
            m_null( i, 0 ) = 0.0;
            m_null( i, 1 ) = 0.0;
            m_null( i, 2 ) = 0.0;
            m_null( i, 3 ) = 0.0;

            missing[ i ] = 1;
        }
    }
}

factor_matrix::factor_matrix(const arma::mat &cov, size_t n)
    : general_matrix( cov, n, 5, 9 )
{
}

void
factor_matrix::update_matrix(const snp_row &row1, const snp_row &row2, arma::uvec &missing)
{
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 && missing[ i ] == 0 )
        {
            double s11 = row1[ i ] == 1 ? 1.0 : 0.0;
            double s12 = row1[ i ] == 2 ? 1.0 : 0.0;
            double s21 = row2[ i ] == 1 ? 1.0 : 0.0;
            double s22 = row2[ i ] == 2 ? 1.0 : 0.0;

            m_alt( i, 0 ) = s11;
            m_alt( i, 1 ) = s12;
            m_alt( i, 2 ) = s21;
            m_alt( i, 3 ) = s22;
            m_alt( i, 4 ) = s11*s21;
            m_alt( i, 5 ) = s11*s22;
            m_alt( i, 6 ) = s12*s21;
            m_alt( i, 7 ) = s12*s22;

            m_null( i, 0 ) = s11;
            m_null( i, 1 ) = s12;
            m_null( i, 2 ) = s21;
            m_null( i, 3 ) = s22;
        }
        else
        {
            m_alt( i, 0 ) = 0.0;
            m_alt( i, 1 ) = 0.0;
            m_alt( i, 2 ) = 0.0;
            m_alt( i, 3 ) = 0.0;
            m_alt( i, 4 ) = 0.0;
            m_alt( i, 5 ) = 0.0;
            m_alt( i, 6 ) = 0.0;
            m_alt( i, 7 ) = 0.0;
            
            m_null( i, 0 ) = 0.0;
            m_null( i, 1 ) = 0.0;
            m_null( i, 2 ) = 0.0;
            m_null( i, 3 ) = 0.0;
            missing[ i ] = 1;
        }
    }
    
}

noia_matrix::noia_matrix(const arma::mat &cov, size_t n)
    : general_matrix( cov, n, 5, 9 )
{
}

void
noia_matrix::update_matrix(const snp_row &row1, const snp_row &row2, arma::uvec &missing)
{
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 && missing[ i ] == 0 )
        {
            double a1 = row1[ i ] - 1;
            double a2 = row2[ i ] - 1;
            double d1 = row1[ i ] == 1 ? 1.0 : 0.0;
            double d2 = row2[ i ] == 1 ? 1.0 : 0.0;
            double aa = a1 * a2;
            double ad = a1 * d2;
            double da = d1 * a2;
            double dd = d1 * d2;

            m_alt( i, 0 ) = a1;
            m_alt( i, 1 ) = a2;
            m_alt( i, 2 ) = d1;
            m_alt( i, 3 ) = d2;
            m_alt( i, 4 ) = a1 * a2;
            m_alt( i, 5 ) = a1 * d2;
            m_alt( i, 6 ) = d1 * a2;
            m_alt( i, 7 ) = d1 * d2;

            m_null( i, 0 ) = a1;
            m_null( i, 1 ) = a2;
            m_null( i, 2 ) = d1;
            m_null( i, 3 ) = d2;
        }
        else
        {
            m_alt( i, 0 ) = 0.0;
            m_alt( i, 1 ) = 0.0;
            m_alt( i, 2 ) = 0.0;
            m_alt( i, 3 ) = 0.0;
            m_alt( i, 4 ) = 0.0;
            m_alt( i, 5 ) = 0.0;
            m_alt( i, 6 ) = 0.0;
            m_alt( i, 7 ) = 0.0;
            
            m_null( i, 0 ) = 0.0;
            m_null( i, 1 ) = 0.0;
            m_null( i, 2 ) = 0.0;
            m_null( i, 3 ) = 0.0;
            missing[ i ] = 1;
        }
    }
    
}

separate_matrix::separate_matrix(const arma::mat &cov, size_t n, separate_mode_t mode)
    : general_matrix( cov, n, 3, 4 )
{
    if( mode == DOM_DOM )
    {
        m_snp1_threshold = m_snp2_threshold = 1;
    }
    else if( mode == REC_DOM )
    {
        m_snp1_threshold = 2;
        m_snp2_threshold = 1;
    }
    else if( mode == DOM_REC )
    {
        m_snp1_threshold = 1;
        m_snp2_threshold = 2;
    }
    else if( mode == REC_REC )
    {
        m_snp1_threshold = m_snp2_threshold = 2;
    }
}

void
separate_matrix::update_matrix(const snp_row &row1, const snp_row &row2, arma::uvec &missing)
{
    for(int i = 0; i < row1.size( ); i++)
    {
        if( row1[ i ] != 3 && row2[ i ] != 3 && missing[ i ] == 0 )
        {
            double snp1 = ( row1[ i ] >= m_snp1_threshold ) ? 1.0 : 0.0;
            double snp2 = ( row2[ i ] >= m_snp2_threshold ) ? 1.0 : 0.0;

            m_alt( i, 0 ) = snp1;
            m_alt( i, 1 ) = snp2;
            m_alt( i, 2 ) = snp1 * snp2;

            m_null( i, 0 ) = snp1;
            m_null( i, 1 ) = snp2;
        }
        else
        {
            m_alt( i, 0 ) = 0.0;
            m_alt( i, 1 ) = 0.0;
            m_alt( i, 2 ) = 0.0;
            
            m_null( i, 0 ) = 0.0;
            m_null( i, 1 ) = 0.0;
            missing[ i ] = 1;
        }
    }
}

model_matrix *
make_model_matrix(const std::string &type, const arma::mat &cov, size_t n)
{
    if( type == "additive" )
    {
        return new additive_matrix( cov, n );
    }
    if( type == "tukey" )
    {
        return new tukey_matrix( cov, n );
    }
    if( type == "factor" )
    {
        return new factor_matrix( cov, n );
    }
    if( type == "noia" )
    {
        return new noia_matrix( cov, n );
    }
    else
    {
        return NULL;
    }
}
