##
# Cleans the given name so that it is suitable for
# a filename.
#
# @param name A possible filename.
#
# @return A suitable filename.
#
def sanitize_name(name):
    return ''.join( x.lower( ) for x in name if x.isalnum( ) )
