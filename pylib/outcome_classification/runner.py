###################################################################
# This script is to run a outcome classifier of droplet impact
# by using deep nueral network. One can train the model by using
# any types of result of the simulation, including spreading diameter,
# droplet velocity, contact angle, and so on. 
# And then, one can get easily categorize the output of your simulation 
# retults into different outcomes.
####################################################################

import argparse
import os
from outcome_classifier import OutcomeClassifier

parser=argparse.ArgumentParser(description='Program to train outcomes of droplet impact on a solid surface')

parser.add_argument('-d', dest='dirstr_filename', type=str, default='DIRSTRUCT', help='Input file for directory structure')
parser.add_argument('-f', dest='analysis_filename', type=str, default='dropsize.out', help='File name containing features in each directory')
parser.add_argument('-c', dest='data_colidx', type=int, default=2, help='The index of the column containing the feature')
parser.add_argument('-t', dest='target_filename', type=str, default='OUTCOMES', help='File name containing targets')
parser.add_argument('-n', dest='nfiles', type=int, default=5, help='The number of files of features` and targets in each directory')
parser.add_argument('-o', dest='outmodel_filename', type=str,default='./classifier_model', help='File name to save models and parameters')
parser.add_argument('-e', dest='end_time', type=int, default=-1, help='Line number to be considered as the end of feature')
parser.add_argument('-s', dest='start_time', type=int, default=0, help='Line number to be considered as the beginning of feature')
parser.add_argument('-p', dest='epochs', type=int, default=1000, help='Line number to be considered as the beginning of feature')

args=vars(parser.parse_args())

#Definition of file structures#
dirstr_filename=args['dirstr_filename']
#File name containing features#
analysis_filename=args['analysis_filename']
#File name contining target values#
target_filename=args['target_filename']
#The number of files for feature in each directory#
nfiles=args['nfiles']
#The column index of the feature in each file
data_colidx=args['data_colidx']
#The name of directory that model is saved
outmodel_filename=args['outmodel_filename']
#The final time of the feature
end_time=args['end_time']
#The first time of the feature
start_time=args['start_time']
#The number of iteration for training
epoch=args['epochs']
#The number of nodes in each layer
num_nodes=[25, 15, 4]


#Defining the class object
outclass=OutcomeClassifier(dirstr_filename, outmodel_filename)

#############################
# For training
#############################

#Structuring the directories for feature
feature_files=outclass.construct_input_struct(analysis_filename, nfiles)
#Structuring the directories for targets
target_files=outclass.construct_target_struct(target_filename)
#Reading features
features=outclass.construct_features(feature_files, data_colidx, start_time, end_time)
#Reading targets
targets=outclass.construct_targets(target_files)
#Defining the DNN model
outclass.dnn_definition(num_nodes)
#Training the DNN
outclass.train(features, targets, epoch)
#Saving parameters
outclass.save_parameters()
#Plotting loss function
outclass.plot_lossfunction()

"""
#############################
# For prediction
#############################
#Loding model
outclass.load_model()
#Structuring the directories for inputs
prediction_set_files=outclass.construct_input_struct(analysis_filename, nfiles)
#Reading prediction sets
prediction_sets=outclass.construct_features(prediction_set_files, data_colidx, start_time, end_time)
#Making prediction
prediction=outclass.output_predict(prediction_sets)
print(prediction)


"""



