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
        self.index = 1

    ##
    # Writes a pair with the next index.
    #
    def write(self, ):
        self.file.write( "snp{0} snp{1}\n".format( self.index, self.index + 1 ) )
        self.index += 2

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )
