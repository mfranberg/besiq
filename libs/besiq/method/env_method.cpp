#include <plink/plink_file.hpp>

#include <besiq/method/env_method.hpp>

void
run_env_method(method_env_type &method, genotype_matrix_ptr genotype_matrix)
{
    std::cout.precision( 4 );
    std::cout << "snp\t";
    method.init( std::cout );
    std::cout << "\tN" << std::endl;

    for(int i = 0; i < genotype_matrix->size( ); i++)
    {
        const snp_row &row = genotype_matrix->get_row( i );
        std::string name = genotype_matrix->get_snp_names( )[ i ];

        std::cout << name << "\t";
        method.run( row, std::cout );
        std::cout << "\t" << method.num_usable_samples( row ) << std::endl;
    }
}
