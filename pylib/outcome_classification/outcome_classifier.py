#######################################################
##   Definition of the class "Outcome Classifier"    ##
#######################################################

import numpy as np
import os
import itertools
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import load_model
from tensorflow.keras.losses import SparseCategoricalCrossentropy
from tensorflow.keras.optimizers import Adam
import matplotlib.pyplot as plt



class OutcomeClassifier:
    #Constructor of the class
    def __init__(self, dirstr_fname, model_fname):
        self.dirstr_fname=dirstr_fname
        self.model_fname=model_fname
        self.outdic={"R":0, "D":1, "N":2, "S":3}
        self.inverse_outdic={val:key for key, val in self.outdic.items()}
    #Structuring directories to be analyzed
    def construct_input_struct(self, analysis_fname, nfiles):
        self.analysis_fname=analysis_fname
        self.nfiles=nfiles
        analfn_pre=analysis_fname[:analysis_fname.find('.')]
        analfn_aft=analysis_fname[analysis_fname.find('.'):]
        try:
            dirstr_file=open(self.dirstr_fname, "r")
        except:
            print(f'[Error] File \'{self.dirstr_fname}\' does not exists.')
        dirstr_lines=dirstr_file.readlines()
        dirstr=[]
        for dirstr_line in dirstr_lines:
            dirstr.append(dirstr_line.split())
        currdir=os.getcwd()
        dir_iter=[list(idx) for idx in itertools.product(*dirstr)]
        allpaths=[]
        for i in range(len(dir_iter)):
            pathname=currdir
            for subdir in dir_iter[i]:
                pathname=pathname+'/'+subdir
            allpaths.append(pathname)

        allfiles=[]
        for path in allpaths:
            for i in range(self.nfiles):
                allfiles.append(path+'/'+analfn_pre+str(i+1)+analfn_aft)


        return allfiles
    #Structuring directories to get targets
    def construct_target_struct(self, target_fname):
        try:
            dirstr_file=open(self.dirstr_fname, "r")
        except:
            print(f'[Error] File \'{self.dirstr_fname}\' does not exists.')
        dirstr_lines=dirstr_file.readlines()
        dirstr=[]
        for dirstr_line in dirstr_lines:
            dirstr.append(dirstr_line.split())
        currdir=os.getcwd()
        dir_iter=[list(idx) for idx in itertools.product(*dirstr)]
        alltargets=[]
        for i in range(len(dir_iter)):
            pathname=currdir
            for subdir in dir_iter[i]:
                pathname=pathname+'/'+subdir
            alltargets.append(pathname+'/'+target_fname)

        return alltargets

    #Construction of features 
    def construct_features(self, feature_files, data_colidx, start_time, end_time):
        x=[]
        for file in feature_files:
            try:
                data=np.loadtxt(file, comments='#')
            except:
                print(f"[Error] File \'{file}\' does not exists.")
            
            x.append(data[start_time:end_time, data_colidx-1])
        x=np.array(x)    
        x=x/x.shape[0]      #normalization
                
        return x

    #Construction of targets
    def construct_targets(self, target_files):
        y=[]
        for file in target_files:
            try:
                data=open(file, "r").readlines()
            except:
                print(f"[Error] File \'{file}\' does not exists.")
            for out in data[:self.nfiles]:
                y.append(self.outdic[out[0]])
        return np.array(y)


    #Defining DNN model
    def dnn_definition(self, num_nodes):
        num_layers=len(num_nodes)
        activation_function=['relu' if i!=num_layers-1 else 'linear' for i in range(num_layers) ]
        self.model=Sequential([Dense(units=num_nodes[i], activation=activation_function[i]) for i in range(num_layers)])
        self.model.compile(loss=SparseCategoricalCrossentropy(from_logits=True),  #<-- Note
            optimizer=Adam(0.001),)
    #Training
    def train(self, features, targets, epochs):
        self.history=self.model.fit(features, targets, epochs=epochs)
        return self.history
    #Saving parameters
    def save_parameters(self):
        self.model.save(self.model_fname)
    #Load model and parameters
    def load_model(self):
        self.model=load_model(self.model_fname)
        self.model.summary()
    #Prediction of the output
    def output_predict(self, inputs):
        pre_predicted=self.model.predict(inputs)
        predicted=[self.inverse_outdic[np.argmax(prob)] for prob in pre_predicted]
        return predicted
    #Plotting loss function
    def plot_lossfunction(self):
        plt.plot(self.history.history['loss'])
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train', 'test'], loc='upper left')
        plt.show()



