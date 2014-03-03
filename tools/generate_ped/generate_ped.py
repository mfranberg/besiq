import argparse

from plink import generate

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser( description='Generates interaction pairs from a given model.' )
    arg_parser.add_argument( '--model', metavar='model', nargs=9, type=generate.probability, help='Space-separated list of floating point numbers that represents the penetrance matrix, specified row-wise from left to right.', required = True )
    arg_parser.add_argument( '--maf', metavar='maf', nargs=2, type=generate.probability, help='Minor allele frequency of the two snps.', required = True )
    arg_parser.add_argument( '--sample-maf',  action="store_true", help='The --maf is treated as a range and maf is sampled uniformly in this range.', default = False )
    arg_parser.add_argument( '--ncases', metavar='ncases', type=int, help='Number of cases.', required = True )
    arg_parser.add_argument( '--ncontrols', metavar='ncontrols', type=int, help='Number of controls.', required = True )
    arg_parser.add_argument( '--npairs', metavar='npairs', type=int, help='Number of interaction pairs', required = True )
    arg_parser.add_argument( '--ld', metavar='ld', type=generate.probability, help='Strength of LD (ignores second maf).' )
    
    arg_parser.add_argument( '--out', metavar='output_file', help='Output .tped file.', required = True )

    args = arg_parser.parse_args( )

    fixed_params = generate.FixedParams( args.maf, args.ld, args.ncases, args.ncontrols, args.sample_maf )
    models = [ ( args.npairs, args.model, False ) ]
    generate.write_data( fixed_params, models, args.out )
