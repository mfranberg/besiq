import subprocess

##
# Generates data for a model without covariates.
#
# @param params General experiment parameters.
# @param model_params Model parameters.
# @param ld Amount of linkage disequilibrium.
# @param output_prefix Path prefix of the plink file to generate.
#
def generate_data(params, model_params, ld, output_prefix):
    script_path = "../tools/generate_ped/generate_ped.py"
    cmd = [ "python", script_path,
            "--maf", str( params.maf[ 0 ] ), str( params.maf[ 1 ] ),
            "--ncases", str( params.sample_size[ 1 ] ),
            "--ncontrols", str( params.sample_size[ 0 ] ),
            "--npairs", str( params.num_pairs ),
            "--out", output_prefix ]
    
    cmd.append( "--model" )
    cmd.extend( list( map( str, model_params ) ) )
    if ld > 0.0:
        cmd.extend( [ "--ld", str( ld ) ] )

    subprocess.call( cmd )

##
# Generates data for a model with covariates that
# are sufficient for generating the disease.
#
# @param params General experiment parameters.
# @param model_params Model parameters.
# @param cov Covariate parameters.
# @param output_prefix Path prefix of the plink file to generate.
#
def generate_sufficient_cov_data(params, model_params, cov, output_prefix):
    script_path = "../tools/generate_ped/generate_sufficient_covariate_ped.py"
    cmd = [ "python", script_path,
            "--maf", str( params.maf[ 0 ] ), str( params.maf[ 1 ] ),
            "--ncases", str( params.sample_size[ 1 ] ),
            "--ncontrols", str( params.sample_size[ 0 ] ),
            "--out", output_prefix ]
    
    cmd.append( "--model" )
    cmd.extend( list( map( str, model_params ) ) )

    if cov:
        cmd.append( "--covariate" )
        cmd.extend( list( map( str, cov ) ) )

    subprocess.call( cmd )


##
# Generates data for a model with covariates from the logistic
# regression model.
#
# @param params General experiment parameters.
# @param snpbeta List of betas for snp1, snp2 and snp1 x snp2.
# @param cov List of covariate parameters as specified to
#            generate_covariate_ped.py.
# @param output_prefix Path prefix of the plink file to generate.
#
def generate_cov_data(params, snpbeta, cov, output_prefix):
    script_path = "../tools/generate_ped/generate_covariate_ped.py"
    cmd = [ "python", script_path,
            "--maf", str( params.maf[ 0 ] ), str( params.maf[ 1 ] ),
            "--ncases", str( params.sample_size[ 1 ] ),
            "--ncontrols", str( params.sample_size[ 0 ] ),
            "--out", output_prefix ]

    cmd.append( "--snpbeta" )
    cmd.extend( list( map( str, snpbeta ) ) )

    if cov:
        cmd.append( "--covariate" )
        cmd.extend( list( map( str, cov ) ) )

    subprocess.call( cmd )

