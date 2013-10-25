##
# Helper class that writes a list of pairs
# and whether they are an interaction or not.
#
class CaseFile:
    ##
    # Constructor.
    #
    # @param path Path to the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )
        self.index = 1

    ##
    # Writes a pair with the next indices and info on
    # whether they are an interaction or not.
    #
    # @param snp1 Index of snp1.
    # @param snp2 Index of snp2.
    # @param is_case Boolean that indicates whether this is a real interaction
    #                or not.
    #
    def write(self,  is_case):
        self.file.write( "snp{0} snp{1} {2}\n".format( self.index, self.index + 1, int( is_case ) ) )
        self.index += 2

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )
