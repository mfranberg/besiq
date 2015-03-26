#include <bayesic/method/method.hpp>
#include <bayesic/pairfile.hpp>

void run_method(method_type &method, const std::vector<snp_row> &genotype_matrix, const std::vector<std::string> &loci, pairfile &pairs)
{
    std::cout.precision( 4 );
    std::cout << "snp1 snp2\t";
    method.init( std::cout );
    std::cout << "\tN" << std::endl;
    
    /* Create a map from snp names in the file to the genotype file */
    std::map<std::string, size_t> snp_to_index;
    for(int i = 0; i < loci.size( ); i++)
    {
        snp_to_index[ loci[ i ] ] = i;
    }

    std::pair<std::string, std::string> pair;
    while( pairs.read( pair ) )
    {
        if( snp_to_index.count( pair.first ) <= 0 || snp_to_index.count( pair.second ) <= 0 )
        {
            continue;
        }

        size_t snp1_index = snp_to_index[ pair.first ];
        size_t snp2_index = snp_to_index[ pair.second ];

        const snp_row &row1 = genotype_matrix[ snp1_index ];
        const snp_row &row2 = genotype_matrix[ snp2_index ];

        std::cout << pair.first << " " << pair.second << "\t";
        method.run( row1, row2, std::cout );
        std::cout << "\t" << method.num_usable_samples( row1, row2 ) << std::endl;
    }
}
