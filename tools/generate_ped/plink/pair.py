##
# Helper class that writes a list of snp pairs.
#
class PairFile:
    ##
    # Constructor.
    #
    # @param path Path to the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )

    ##
    # Writes a pair with the given indicies ot
    # the file.
    #
    # @param snp1 Index of snp1.
    # @param snp2 Index of snp2.
    #
    def write(self, snp1, snp2):
        self.file.write( "snp{0} snp{1}\n".format( snp1, snp2 ) )

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )
