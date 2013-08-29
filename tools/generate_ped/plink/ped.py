##
# Simple helper class for sequentially writing .ped files.
#
#
class PedFile:
    ##
    # Constructor.
    #
    # @param path Path to the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )
        self.fam_id = 1

    ##
    # Write an individual, the iid and fid will be incremented for
    # each write. 
    #
    # @param alleles A list representing the allele at each loci.
    # @param phenotype The phenotype of the individual.
    #
    def write(self, alleles, phenotype):
        self.file.write( "{0} {0} 0 0 1 {1} {2}\n".format( self.fam_id, phenotype, " ".join( alleles ) ) )
        self.fam_id += 1

    ##
    # Close the file.
    #
    def close(self):
        self.file.close( )

