##
# Helper class that writes covariates to a file.
#
class CovFile:
    ##
    # Constructor.
    #
    # @param path Path to the file.
    # @param covariate_names Names of the covariates.
    #
    def __init__(self, path, covariate_names):
        self.file = open( path, "w" )
        self.file.write( "{0}\n".format( "\t".join( covariate_names ) ) )
        self.iid = 0

    ##
    # Writes the covariates for and individual.
    #
    # @param covariate_values Values of the covariates.
    #
    def write(self, covariate_values):
        self.file.write( "{0}\n".format( "\t".join( map( str, covariate_values ) ) ) )
        self.iid += 1

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )
