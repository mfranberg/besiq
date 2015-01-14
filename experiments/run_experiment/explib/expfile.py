import os
from . import util

##
# Handles active output files for different methods.
#
class MethodHandler:
    ##
    # @param base_dir Base output directory.
    # 
    def __init__(self, base_dir):
        ##
        # Base output directory.
        #
        self.base_dir = base_dir

        ##
        # Name of the experiment.
        #
        self.experiment_name = ""

        ##
        # Id of the model.
        #
        self.model_id = ""

        ##
        # Additional postfix for method output filenames.
        #
        self.postfix = ""

        ##
        # Mapping from method to output file.
        #
        self.method_files = dict( )

    ##
    # Resets the object so that it can handle a new set of methods.
    #
    # @param experiment_name Name of the experiment.
    # @param model_id Id of the model.
    # @param postfix Additional postfix for method output filenames.
    #
    def start_experiment(self, experiment_name, model_id, postfix = ""):
        self.experiment_name = experiment_name
        self.model_id = model_id
        self.postfix = "_" + postfix
        self.method_files = dict( )
 
    ##
    # Seeks to the beginning of all files, use when output is finished
    # and you want to read the files.
    #
    def reset_files(self):
        for method_file in self.method_files.itervalues( ):
            method_file.seek( 0 )

    ##
    # Generates an output path for the method with the given name.
    #
    # @param method_name Name of the method.
    #
    # @return Returns a generated method output path.
    #
    def generate_output_path(self, method_name):
        corrected_method_name = util.sanitize_name( method_name )
        corrected_experiment_name = util.sanitize_name( self.experiment_name )
        output_file = "{0}_{1}_effect{2}{3}.out".format( corrected_experiment_name,
                                                      corrected_method_name,
                                                      self.model_id,
                                                      self.postfix )
        return os.path.join( self.base_dir, output_file )

    ##
    # Returns the output file corresponding to the given method name.
    #
    # @param method_name Name of the method.
    #
    # @return The output file corresponding to the given method name.
    #
    def get_output_file(self, method_name):
        method_file = self.method_files.get( method_name, None )
        if method_file:
            return method_file
        else:
            output_path = self.generate_output_path( method_name )
            self.method_files[ method_name ] = open( output_path, "w+" )
            return self.method_files[ method_name ]

##
# Handles various paths in the experiment directory.
#
class PathHandler:
    ##
    # @param output_dir Path to the output directory.
    # 
    def __init__(self, output_dir, subdir_list):
        ##
        # Experiment output directory.
        #
        self.output_dir = output_dir

        ##
        # List of subdirectories to create.
        #
        self.subdirs = dict( )
        for subdir in subdir_list:
            self.subdirs[ subdir ] = os.path.join( self.output_dir, subdir )

    ##
    # Returns the list of subdirs this object handles.
    #
    # @return The list of subdirs this object handles.
    #
    def get_subdirs(self):
        return self.subdirs.values( )
    
    ##
    # Returns the full path to the generated plink file.
    #
    # @return The full path to the generated plink file.
    #
    def get_plink_path(self):
        return os.path.join( self.subdirs[ "tmp" ], "plink" )

    ##
    # Returns the full path to the method output subdirectory.
    #
    # @return The full path to the method output subdirectory.
    #
    def get_method_base_path(self):
        return self.subdirs[ "method" ]
 
    ##
    # Returns the full path to the file that contains power estimates.
    #
    # @return The full path to the file that contains power estimates.
    #
    def get_power_path(self, experiment_name):
        corrected_experiment_name = util.sanitize_name( experiment_name )

        return os.path.join( self.subdirs[ "power" ], "power_{0}.out".format( corrected_experiment_name ) )

    ##
    # Returns the full path to the plot file for the experiment.
    #
    # @return The full path to the plot file for the experiment.
    #
    def get_plot_path(self, experiment_name):
        corrected_experiment_name = util.sanitize_name( experiment_name )
        return os.path.join( self.subdirs[ "plots" ], "{0}.pdf".format( corrected_experiment_name ) )
