##
# Helper class that writes a list of pairs
# and the index of the model that they are
# generated from.
#
class ModelFile:
    ##
    # Constructor.
    #
    # @param path Path to the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )
        self.index = 1

    ##
    # Writes a pair with the next indices and the index
    # of the model.
    #
    # @param snp1 Index of snp1.
    # @param snp2 Index of snp2.
    # @param model_index Index of the model.
    #
    def write(self, model_index):
        self.file.write( "snp{0} snp{1} {2}\n".format( self.index, self.index + 1, int( model_index ) ) )
        self.index += 2

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )
