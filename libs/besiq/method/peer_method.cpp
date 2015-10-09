#include <dcdflib/libdcdf.hpp>
#include <besiq/method/peer_method.hpp>
#include <besiq/stats/snp_count.hpp>

peer_method::peer_method(method_data_ptr data)
: method_type::method_type( data )
{
    m_weight = arma::ones<arma::vec>( data->phenotype.size( ) );
}

std::vector<std::string>
peer_method::init()
{
    std::vector<std::string> header;
    header.push_back( "Z_dd_ld" );
    header.push_back( "P_dd_con" );
    header.push_back( "Z_rd_ld" );
    header.push_back( "P_rd_con" );
    header.push_back( "Z_dr_ld" );
    header.push_back( "P_dr_con" );
    header.push_back( "Z_rr_ld" );
    header.push_back( "P_rr_con" );

    return header;
}

/**
 * rr | rd
 * dr | dd
 *
 * 00    | 00/01 | 01
 * 10/00 | 00/11 | 01/11
 * 10    | 10/11 | 11
 *
 * @param counts The counts for each genotype.
 * @param snp1_onestart The minimum number of minor alleles that should be considered a 1,
 *                      at locus 1.
 * @param snp2_onestart The minimum number of minor alleles that should be considered a 1.
 *                      at locus 2.
 */
arma::mat
peer_method::encode_counts(const arma::mat &counts, size_t snp1_onestart, size_t snp2_onestart)
{
    int dd = ( snp1_onestart == 1 ) && ( snp2_onestart == 1 );
    int rd = ( snp1_onestart == 1 ) && ( snp2_onestart == 2 );
    int dr = ( snp1_onestart == 2 ) && ( snp2_onestart == 1 );
    int rr = ( snp1_onestart == 2 ) && ( snp2_onestart == 2 );

    arma::mat encoded_counts( 4, 2 );
    for(int i = 0; i <= 1; i++)
    {
        encoded_counts( 0, i ) = counts( 0, i ) + (rd + rr) * counts( 1, i ) + (dr + rr) * counts( 3, i ) + rr * counts( 4, i );
        encoded_counts( 1, i ) = (dr + dd) * counts( 1, i ) + counts( 2, i ) + dr * counts( 4, i ) + (dr + rr) * counts( 5, i );
        encoded_counts( 2, i ) = (rd + dd) * counts( 3, i ) + rd * counts( 4, i ) + counts( 6, i ) + (rd + rr) * counts( 7, i );
        encoded_counts( 3, i ) = dd * counts( 4, i ) + (rd + dd) * counts( 5, i ) + (dr + dd) * counts( 7, i ) + counts( 8, i );
    }

    return encoded_counts;
}

void
peer_method::compute_ld_p(const arma::mat &counts, float *ld_case_z, float *ld_contrast_z)
{
    double N = arma::accu( counts );
    double N_controls = sum( counts.col( 0 ) );
    double N_cases = sum( counts.col( 1 ) );

    double p_a_case = ( counts( 1, 1 ) + counts( 3, 1 ) ) / N_cases;
    double p_a_control = ( counts( 1, 0 ) + counts( 3, 0 ) ) / N_controls;
    double p_b_case = ( counts( 2, 1 ) + counts( 3, 1 ) ) / N_cases;
    double p_b_control = ( counts( 2, 0 ) + counts( 3, 0 ) ) / N_controls;
    
    double p_ab_case = counts( 3, 1 ) / N_cases;
    double p_ab_control = counts( 3, 0 ) / N_controls;

    double delta_case = ( p_ab_case - p_a_case * p_b_case );
    double delta_control = ( p_ab_control - p_a_control * p_b_control );

    double sigma2_case = ( p_a_case * (1 - p_a_case) * p_b_case * (1 - p_b_case ) ) / N_cases; 
    double sigma2_control = ( p_a_control * (1 - p_a_control) * p_b_control * (1 - p_b_control ) ) / N_controls;
    double sigma_diff = sqrt( sigma2_case + sigma2_control );

    *ld_case_z = delta_case / sqrt( sigma2_case );
    *ld_contrast_z = ( delta_case - delta_control ) / sigma_diff;
}

double
peer_method::run(const snp_row &row1, const snp_row &row2, float *output)
{
    arma::mat counts = joint_count( row1, row2, get_data( )->phenotype, m_weight );
    
    if( arma::min( arma::min( counts ) ) < METHOD_SMALLEST_CELL_SIZE_BINOMIAL )
    {
        return -9;
    }
    set_num_ok_samples( (size_t) arma::accu( counts ) );

    size_t output_index = 0;
    for(int i = 1; i <= 2; i++)
    {
        for(int j = 1; j <= 2; j++)
        {
            arma::mat encoded_counts = encode_counts( counts, i, j );

            float ld_case_z;
            float ld_contrast_z;
            compute_ld_p( encoded_counts, &ld_case_z, &ld_contrast_z );

            output[ output_index++ ] = ld_case_z;
            output[ output_index++ ] = 1 - norm_cdf( ld_contrast_z, 0.0, 1.0 );
        }
    }

    return min_na( min_na( output[ 1 ], output[ 3 ] ), min_na( output[ 5 ], output[ 7 ] ) );
}
