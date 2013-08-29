##
# Helper class for writing .map files.
#
class MapFile:
    ##
    # Constructor.
    #
    # @param path Path for the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )
        self.snp_id = 0

    ##
    # Writes a snp where the name and position is incremented
    # for each write.
    #
    def write(self):
        self.file.write( "1 snp{0} 0 {0}\n".format( self.snp_id ) )
        self.snp_id += 1

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )

