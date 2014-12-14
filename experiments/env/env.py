import argparse
import subprocess
import os

def generate_plink(plink_path, maf, sample_size, num_variants):
    wgas_path = plink_path + ".wgas"
    wgas_file = open( wgas_path, "w" )
    wgas_file.write( "{1} null {0} {0} 1.00 1.00\n".format( maf, num_variants ) )
    wgas_file.close( )

    cmd = [ "plink2",
            "--simulate", wgas_path,
            "--make-bed",
            "--simulate-ncases", str( sample_size ),
            "--simulate-ncontrols", "0",
            "--out", plink_path ]

    print " ".join( cmd )
    subprocess.check_call( cmd, stdout = open( os.devnull, "w" ) )


def generate_pheno(plink_path, model, std, env_path, pheno_path):
    script_path = os.path.dirname( os.path.realpath( __file__ ) )
    generate_path = os.path.join( script_path, "../../tools/generate_ped/generate_env.py" )

    cmd = [ "python",
           generate_path,
           "--env", "0.5", "0.5",
           "--env-file", env_path,
           "--pheno-file", pheno_path,
           "--std", str( std ),
           plink_path ]

    cmd.append( "--model" )
    cmd.extend( map( str, model ) )

    print " ".join( cmd )
    subprocess.check_call( cmd )

def run_test(plink_path, env_path, pheno_path, result_path):
    cmd = [ "bayesic-env",
            "-m", "stepwise",
            "-p", pheno_path,
            env_path,
            plink_path ]
    with open( result_path, "w" ) as result_file:
        print " ".join( cmd )
        subprocess.check_call( cmd, stdout = result_file )

def count_lines(result_path, level):
    num_lines = 0
    for line in open( result_path, "r" ):
        try:
            p = float( column[ 2 + level ] )
            num_lines = num_lines + 1
        except:
            continue

    return num_lines

def update_final(result_path, final_file):
    prev_num_tests = count_lines( result_path, 0 )
    for level in range( 4 ):
        prev_path = result_path + ".level" + str( level - 1 )
        if level == 0:
            prev_path = result_path

        level_file = open( result_path + ".level" + str( level ), "w" )
        cur_num_tests = 0
        for line in open( prev_path, "r" ):
            column = line.strip( ).split( )
            try:
                p = float( column[ 2 + level ] )
                if p < 0.05 / prev_num_tests:
                    level_file.write( line )
                cur_num_tests += 1
            except:
                continue

        level_file.close( )
        prev_num_tests = cur_num_tests

    final_file.write( "{0}\n".format( prev_num_tests ) )

def main():
    parser = argparse.ArgumentParser( description = "Runs bayesic-env repeatedly and calculates the number of successful inferences." )
    parser.add_argument( 'output_dir', type = str, help = 'The output directory.' )
    parser.add_argument( '--plink-file', type = str, help = 'Plink prefix.', default = None )
    parser.add_argument( '-n', type = int, help = 'Number of replicates.', default = 100 )
    parser.add_argument( '--model', type = float, nargs = 6, help = 'Gene environment model (as passed to generate_env.py)', required = True )
    parser.add_argument( '--std', type = float,  help = 'Standard deviation.', required = True )
    parser.add_argument( '--maf', type = float,  help = 'Minor allele frequency if plink file not set', default = 0.3 )
    parser.add_argument( '--num-variants', type = int,  help = 'Number of variants.', default = 10000 )
    parser.add_argument( '--sample-size', type = int,  help = 'Number of samples if plink file not set', default = 3000 )

    args = parser.parse_args( )

    plink_path = args.plink_file
    if not plink_path:
        plink_path = os.path.join( args.output_dir, "plink" )
        generate_plink( plink_path, args.maf, args.sample_size, args.num_variants )

    if not os.path.isdir( args.output_dir ):
        os.makedirs( args.output_dir )

    final_file = open( os.path.join( args.output_dir, "final.out" ), "w" )
    for i in range( args.n ):
        pheno_path = os.path.join( args.output_dir, "plink.pheno" )
        env_path = os.path.join( args.output_dir, "plink.env" )
        generate_pheno( plink_path, args.model, args.std, env_path, pheno_path )

        result_path = os.path.join( args.output_dir, "bayesic.out" )
        run_test( plink_path, env_path, pheno_path, result_path )

        update_final( result_path, final_file )

    final_file.close( )

if __name__ == "__main__":
    main( )
