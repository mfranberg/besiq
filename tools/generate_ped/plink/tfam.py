##
# Helper class for writing .tfam files.
#
class TFamFile:
    ##
    # Constructor.
    #
    # @param path Path for the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )
        self.iid = 1

    ##
    # Writes a sample with the given phenotypes. FID
    # and IID is incremented for each write.
    #
    # @param phenotype Phenotype of the individual
    #
    def write(self, phenotype):
        self.file.write( "{0} {0} 0 0 1 {1}\n".format( self.iid, phenotype ) )
        self.iid += 1

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )

