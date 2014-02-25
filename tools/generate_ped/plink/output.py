from .tfam import TFamFile
from .tped import TPedFile
from .pair import PairFile
from .case import CaseFile
from .model import ModelFile

class OutputFiles:
    def __init__(self, path):
        self.tped_file = TPedFile( path + ".tped" )
        self.tfam_file = TFamFile( path + ".tfam" )
        self.pair_file = PairFile( path + ".pair" )
        self.case_file = CaseFile( path + ".case" )
        self.model_file = ModelFile( path + ".model" )

    def close(self):
        self.tped_file.close( )
        self.tfam_file.close( )
        self.pair_file.close( )
        self.case_file.close( )
        self.model_file.close( )
