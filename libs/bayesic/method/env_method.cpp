#include <bayesic/method/env_method.hpp>

void
run_env_method(method_env_type &method, const std::vector<snp_row> &genotype_matrix, const std::vector<std::string> &loci)
{
    std::cout.precision( 4 );
    std::cout << "snp\t";
    method.init( std::cout );
    std::cout << "\tN" << std::endl;

    for(int i = 0; i < loci.size( ); i++)
    {
        const snp_row &row = genotype_matrix[ i ];
        std::string name = loci[ i ];

        std::cout << name << "\t";
        method.run( row, std::cout );
        std::cout << "\t" << method.num_usable_samples( row ) << std::endl;
    }
}
