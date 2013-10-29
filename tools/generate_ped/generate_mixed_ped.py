import argparse

from plink import generate

##
# Parses a set of models in a given model file.
#
# @param model_file A file where each row contains: number of pairs,
#                   boolean, and 9 penetrance values.
#
# @return An iterator over the models in the file.
#
def parse_models(model_file):
    for line in model_file:
        columns = line.strip( ).split( )

        if len( columns ) != 11:
            print( "Could not parse line in model file." )
            exit( 1 )

        try:
            num_pairs = int( columns[ 0 ] )
            is_case = bool( int( columns[ 1 ], 2 ) )
            params = list( map( float, columns[ 2: ] ) )

            yield (num_pairs, params, is_case)

        except ValueError:
            print( "Could not parse value in model file." )
            exit( 1 )

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser( description='Generates interaction pairs from a given model.' )
    arg_parser.add_argument( '--model', metavar='model', type=argparse.FileType( "r" ), help='File where each row contains the number of pairs and 9 penetrance values.', required = True )
    arg_parser.add_argument( '--maf', metavar='maf', nargs=2, type=generate.probability, help='Minor allele frequency of the two snps.', required = True )
    arg_parser.add_argument( '--ncases', metavar='ncases', type=int, help='Number of cases.', required = True )
    arg_parser.add_argument( '--ncontrols', metavar='ncontrols', type=int, help='Number of controls.', required = True )
    arg_parser.add_argument( '--ld', metavar='ld', type=generate.probability, help='Strength of LD (ignores second maf).' )
    
    arg_parser.add_argument( '--out', metavar='output_file', help='Output .tped file.', required = True )

    args = arg_parser.parse_args( )

    fixed_params = generate.FixedParams( args.maf, args.ld, args.ncases, args.ncontrols )
    models = parse_models( args.model )
    generate.write_data( fixed_params, models, args.out )
