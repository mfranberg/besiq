##
# Simple helper class for sequentially writing .tped files
# with respect to the loci.
#
#
class TPedFile:
    ##
    # Constructor.
    #
    # @param path Path to the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )
        self.snp_id = 0

    ##
    # Write a loci, the name of the loci will be snp followed
    # by an incremented number.
    #
    # @param alleles A list representing the allele for each individual
    #
    def write(self, alleles):
        self.file.write( "1 snp{0} 0 {0} {1}\n".format( self.snp_id, " ".join( alleles ) ) )
        self.snp_id += 1

    ##
    # Close the file.
    #
    def close(self):
        self.file.close( )

