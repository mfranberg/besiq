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

    if params.sample_maf:
        cmd.append( "--sample-maf" )
    
    cmd.append( "--model" )
    cmd.extend( list( map( str, model_params ) ) )
    if ld > 0.0:
        cmd.extend( [ "--ld", str( ld ) ] )

    subprocess.call( cmd )

##
# Generates data for a mix of models.
#
# @param params General experiment parameters.
# @param model_params A list of tuples containing: the number
#                     of models to generate, whether it is an interaction,
#                     the parameters of the model as in generate_data.
# @param output_prefix Path prefix of the plink files to generate.
#
def generate_mixed_data(params, model_params, output_prefix):
    model_file_path = output_prefix + ".models"
    with open( model_file_path, "w" ) as model_file:
        for num_tests, is_case, m_params in model_params:
            param_str = " ".join( map( str, m_params ) )
            model_file.write( "{0} {1} {2}\n".format( num_tests, is_case, param_str ) )
   
    script_path = "../tools/generate_ped/generate_mixed_ped.py"
    cmd = [ "python", script_path,
            "--maf", str( params.maf[ 0 ] ), str( params.maf[ 1 ] ),
            "--ncases", str( params.sample_size[ 1 ] ),
            "--ncontrols", str( params.sample_size[ 0 ] ),
            "--model", model_file_path,
            "--out", output_prefix ]

    subprocess.call( cmd )

##
# Generates data for all interaction models.
#
# @param params General experiment parameters.
# @param h Heritability.
# @param base_risk The population risk.
# @param output_prefix Path prefix of the plink files to generate.
#
def generate_all_data(params, h, base_risk, output_prefix): 
    script_path = "../tools/generate_ped/generate_all.py"
    cmd = [ "python", script_path,
            "--maf", str( params.maf[ 0 ] ), str( params.maf[ 1 ] ),
            "--ncases", str( params.sample_size[ 1 ] ),
            "--ncontrols", str( params.sample_size[ 0 ] ),
            "--base-risk", str( base_risk ),
            "--heritability", str( h ),
            "--num-pairs", str( params.num_pairs ),
            "--out", output_prefix ]

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

