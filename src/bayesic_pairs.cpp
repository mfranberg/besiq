#include <iostream>
#include <string>
#include <vector>

#include <gzstream/gzstream.hpp>

#include <plink/plink_file.hpp>
#include <cpp-argparse/OptionParser.h>

#include <bayesic/stats/snp_count.hpp>

using namespace optparse;

const std::string USAGE = "bayesic-pairs genotype_plink_prefix";
const std::string VERSION = "bayesic-pairs 1.0.0";
const std::string DESCRIPTION = "Generates a list of interactions to test with Bayesic.";
const std::string EPILOG = "";

/**
 * Computes the minor allele frequency for each
 * snp in the given plink file.
 *
 * @param genotype_file Plink file.
 *
 * @return A vector containing the maf of all snps.
 */
std::vector<double>
compute_maf(plink_file_ptr &genotype_file)
{
    std::vector<double> maf_vec;
    snp_row row;
    while( genotype_file->next_row( row ) )
    {
        double maf = compute_real_maf( row );
        if( maf > 0.5 )
        {
            maf = 1.0 - maf;
        }

        maf_vec.push_back( maf );
    }

    return maf_vec;
}

/**
 * A vector of pairs that represents pairs of genes.
 */
typedef std::vector< std::pair<std::string, std::string> > pair_vector;

/**
 * Parses a file in which each line contains two names.
 *
 * @param path Path of the file.
 *
 * @return The parsed gene names.
 */
pair_vector
parse_genes(const std::string &path)
{
    pair_vector pairs;
    std::ifstream gene_file( path.c_str( ) );
    while( gene_file.good( ) )
    {
        std::string gene1;
        std::string gene2;

        if( ( gene_file >> gene1 ) && ( gene_file >> gene2 ) )
        {
            pairs.push_back( std::make_pair( gene1, gene2 ) );
        }
    }

    return pairs;
}

/**
 * Creates an opposite map of a vector, which is 
 * indexed by the value and maps to the key, in this
 * case the index of the vector.
 *
 * @param loci A list of locus names.
 *
 * @return A map from locus name to its index.
 */
std::map< std::string, size_t >
create_loci_index(const std::vector<std::string> &loci)
{
    std::map< std::string, size_t > index;
    for(int i = 0; i < loci.size( ); i++)
    {
        index[ loci[ i ] ] = i;
    }

    return index;
}

/**
 * Parses a file that contains in which each line contains a
 * gene name and a locus. Effectively grouping each snp in genes.
 *
 * @param path Path to the file.
 * @param loci A list of locus names.
 *
 * @return A map from gene name to a list of loci.
 */
std::map< std::string, std::vector<size_t> >
parse_gene_locus(const std::string &path, const std::vector<std::string> &loci)
{
    std::map< std::string, std::vector<size_t> > gene_locus;
    std::map< std::string, size_t > index = create_loci_index( loci );
    std::ifstream gene_locus_file( path.c_str( ) );
    while( gene_locus_file.good( ) )
    {
        std::string gene;
        std::string locus;

        if( ( gene_locus_file >> gene ) && ( gene_locus_file >> locus ) )
        {
            gene_locus[ gene ].push_back( index[ locus ] );
        }
    }

    return gene_locus;
}

/**
 * A group of options that determine which pairs of snps 
 * that should be outputted.
 */
struct output_options
{
    /**
     * The locus names from the genotype files.
     */
    std::vector<std::string> loci;

    /**
     * The list of maf for each snp.
     */
    std::vector<double> maf_vec;

    /**
     * The threshold that determine whether a snp should
     * be included or not.
     */
    double maf_threshold;

    /**
     * The threshold that determines whether two snps should
     * be included or not, based on the product of their mafs.
     */
    double combined_threshold;
};

/**
 * For each gene, outputs all pairs of snps in that gene.
 *
 * @param output Output stream.
 * @param oo Output options.
 * @param gene_locus Map from gene to the loci belonging to that gene.
 */
void output_within(std::ostream &output, const output_options &oo, const std::map< std::string, std::vector<size_t> > &gene_locus)
{
    std::map< std::string, std::vector<size_t> >::const_iterator it;
    for(it = gene_locus.begin( ); it != gene_locus.end( ); ++it)
    {
        const std::vector<size_t> &indices = it->second;
        for(int i = 0; i < indices.size( ); i++)
        {
            int snp1 = indices[ i ];
            if( oo.maf_vec[ snp1 ] < oo.maf_threshold )
            {
                continue;
            }

            for(int j = i + 1; j < indices.size( ); j++)
            {
                int snp2 = indices[ j ];
                if( oo.maf_vec[ snp2 ] >= oo.maf_threshold && (oo.maf_vec[ snp1 ] * oo.maf_vec[ snp2 ]) >= oo.combined_threshold )
                {
                    output << oo.loci[ snp1 ] << " " << oo.loci[ snp2 ] << "\n";
                }
            }
        }
    }
}

/**
 * For each pair of genes, outputs all pair of snps in those genes.
 *
 * @param output Output stream.
 * @param oo Output options.
 * @param gene_locus Map from gene to the loci belonging to that gene.
 */
void output_between(std::ostream &output, const output_options &oo, std::map< std::string, std::vector<size_t> > &gene_locus)
{
    std::vector<std::string> genes;
    std::map< std::string, std::vector<size_t> >::const_iterator it;
    for(it = gene_locus.begin( ); it != gene_locus.end( ); ++it)
    {
        genes.push_back( it->first );
    }

    for(int g1 = 0; g1 < genes.size( ); g1++)
    {
        const std::vector<size_t> &indices1 = gene_locus[ genes[ g1 ] ];
        for(int g2 = g1 + 1; g2 < genes.size( ); g2++)
        {
            const std::vector<size_t> &indices2 = gene_locus[ genes[ g2 ] ];
            
            for(int i = 0; i < indices1.size( ); i++)
            {
                int snp1 = indices1[ i ];
                if( oo.maf_vec[ snp1 ] < oo.maf_threshold )
                {
                    continue;
                }

                for(int j = 0; j < indices2.size( ); j++)
                {
                    int snp2 = indices2[ j ];
                    if( oo.maf_vec[ snp2 ] >= oo.maf_threshold && (oo.maf_vec[ snp1 ] * oo.maf_vec[ snp2 ]) >= oo.combined_threshold )
                    {
                        output << oo.loci[ snp1 ] << " " << oo.loci[ snp2 ] << "\n";
                    }
                }
            }
        }
    }
}

/**
 * For each pair of genes in the given list, outputs all pair of snps in those genes.
 *
 * @param output Output stream.
 * @param oo Output options.
 * @param gene_locus Map from gene to the loci belonging to that gene.
 * @param gene_gene A list of pairs of genes to be considered.
 */
void output_between_restrict(std::ostream &output, const output_options &oo, std::map< std::string, std::vector<size_t> > &gene_locus, const pair_vector &gene_gene)
{
    for(int g = 0; g < gene_gene.size( ); g++)
    {
        const std::pair<std::string, std::string> &genes = gene_gene[ g ];
        
        const std::vector<size_t> &indices1 = gene_locus[ genes.first ];
        const std::vector<size_t> &indices2 = gene_locus[ genes.second ];
        
        for(int i = 0; i < indices1.size( ); i++)
        {
            int snp1 = indices1[ i ];
            if( oo.maf_vec[ snp1 ] < oo.maf_threshold )
            {
                continue;
            }

            for(int j = 0; j < indices2.size( ); j++)
            {
                int snp2 = indices2[ j ];
                if( oo.maf_vec[ snp2 ] >= oo.maf_threshold && (oo.maf_vec[ snp1 ] * oo.maf_vec[ snp2 ]) >= oo.combined_threshold )
                {
                    output << oo.loci[ snp1 ] << " " << oo.loci[ snp2 ] << "\n";
                }
            }
        }

    }
}

/**
 * Outputs all pairs of snps.
 *
 * @param output Output stream.
 * @param oo Output options.
 */
void output_all(std::ostream &output, const output_options &oo)
{
    for(int i = 0; i < oo.loci.size( ); i++)
    {
        if( oo.maf_vec[ i ] < oo.maf_threshold )
        {
            continue;
        }

        for(int j = i + 1; j < oo.loci.size( ); j++)
        {
            if( oo.maf_vec[ j ] >= oo.maf_threshold && (oo.maf_vec[ i ] * oo.maf_vec[ j ]) >= oo.combined_threshold )
            {
                output << oo.loci[ i ] << " " << oo.loci[ j ] << "\n";
            }
        }
    }
}

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG );
    
    parser.add_option( "-m", "--maf" ).type( "float" ).set_default( 0.0 ).help( "Remove pairs where one of the SNPs have a maf less than this." );
    parser.add_option( "-c", "--combined-maf" ).type( "float" ).set_default( 0.0 ).help( "Remove pairs where the product of the MAFs is less than this." );
    parser.add_option( "-w", "--within" ).help( "Only output pairs of snps within the genes given by this file." );
    parser.add_option( "-b", "--between" ).help( "Only output pairs of snps between pairs of genes as specified by this file." );
    parser.add_option( "-r", "--restrict" ).help( "Used with --between to only check the pair of genes in this list." );
    parser.add_option( "-o", "--out" ).help( "Name of the output file, will be gzipped." ).set_default( "" );

    Values options = parser.parse_args( argc, argv );
    std::vector<std::string> args = parser.args( );
    if( args.size( ) != 1 )
    {
        printf( "bayesic-pairs: error: Pairs or genotypes is missing.\n" );
        parser.print_help( );
        exit( 1 );
    }

    std::ios_base::sync_with_stdio( false );
    gz::ogzstream output_file( ( (std::string) options.get( "out" ) ).c_str( ) );
    std::ostream &output = options.is_set( "out" ) ? output_file : std::cout;

    output_options oo;
    oo.maf_threshold = (double) options.get( "maf" );
    oo.combined_threshold = (double) options.get( "combined_maf" );

    plink_file_ptr genotype_file = open_plink_file( args[ 0 ] );
    oo.maf_vec = compute_maf( genotype_file );
    oo.loci = genotype_file->get_locus_names( );
    if( options.is_set( "within" ) )
    {
        std::map< std::string, std::vector<size_t> > gene_locus = parse_gene_locus( options[ "within" ].c_str( ), oo.loci );
        output_within( output, oo, gene_locus );
    }
    else if( options.is_set( "between" ) )
    {
        std::map< std::string, std::vector<size_t> > gene_locus = parse_gene_locus( options[ "between" ].c_str( ), oo.loci );
        if( !options.is_set( "restrict" ) )
        {
            output_between( output, oo, gene_locus );
        }
        else
        {
            pair_vector gene_pairs = parse_genes( options[ "restrict" ].c_str( ) );
            output_between_restrict( output, oo, gene_locus, gene_pairs );
        }
    }
    else
    {
        output_all( output, oo );
    }

    return 1;
}
