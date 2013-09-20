#include <method/method.hpp>

void run_method(method_type &method, const std::vector<snp_row> &genotype_matrix, const std::vector<std::string> &loci, pair_iter &pairs)
{
    std::cout.precision( 4 );
    std::cout << "snp1 snp2\t";
    method.init( std::cout );
    std::cout << "\tN" << std::endl;

    std::pair<size_t, size_t> pair;
    while( pairs.get_pair( &pair ) )
    {
        const snp_row &row1 = genotype_matrix[ pair.first ];
        const snp_row &row2 = genotype_matrix[ pair.second ];

        std::string name1 = loci[ pair.first ];
        std::string name2 = loci[ pair.second ];

        std::cout << name1 << " " << name2 << "\t";
        method.run( row1, row2, std::cout );
        std::cout << "\t" << method.num_usable_samples( row1, row2 ) << std::endl;
    }
}
