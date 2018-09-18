"""
Created:        12 November  2016
Last Updated:   25 February  2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Script for running the deep learning implementation

To run:
$ python python/runDeepLearning.py config/mlconfig.txt
-- the second argument is the text file with configurations for the NN/setup

 {"none",0},    // :: NONE = QCD (background)
 {"QB",1},      // :: QB-Q = Signal AK8(QB) + AK4(Q)
 {"W",2},       // :: QQ-B = Signal AK8(W)  + AK4(B)
 {"other",3"    // :: QB/QQ+AK4 = Background ("signal-like")
"""
import os
import sys
import json
from collections import Counter
from time import strftime,localtime


import plotlabels as plb
cwd = os.getcwd()
hpd = cwd.rstrip("goldilocks")+"/asimov/python/"
if hpd not in sys.path:
    sys.path.insert(0,hpd)

import util
from training import Training
from config import Config


print
print " ------------------------------ "
print " *  Goldilocks Deep Learning  * "
print " ------------------------------ "
print

date = strftime("%d%b%Y-%H%M")
vb   = util.VERBOSE()

## Set configuration options ##
config = Config(sys.argv[1])
vb.level = config.verbose_level
vb.initialize()

if not config.runTraining and not config.runInference:
    vb.ERROR("RUN :  No configuration set ")
    vb.ERROR("RUN :  Please set the arguments 'runTraining' or 'runInference' to define workflow ")
    vb.ERROR("RUN :  Exiting.")
    sys.exit(1)


## Setup Deep Learning class
dnn = Training()

dnn.variable_labels = plb.variable_labels()
dnn.sample_labels   = plb.sample_labels()

dnn.hep_data   = config.hep_data
dnn.model_name = config.dnn_data
dnn.msg_svc    = vb
dnn.treename   = config.treename
dnn.useLWTNN   = True
dnn.dnn_name   = "dnn"
dnn.output_dim = config.output_dim
dnn.loss       = config.loss
dnn.init       = config.init
dnn.nNodes     = config.nNodes
dnn.dropout    = None
dnn.metrics    = config.metrics
dnn.features   = config.features
dnn.epochs     = config.epochs
dnn.optimizer  = config.optimizer
dnn.input_dim  = len(config.features)
dnn.batch_size = config.batch_size
dnn.activations    = config.activation.split(',')
dnn.nHiddenLayers  = config.nHiddenLayers
dnn.earlystopping  = {'monitor':'loss','min_delta':0.0001,'patience':10,'mode':'auto'}
dnn.runDiagnostics = True
dnn.classes = {"multijet":0,"BQ":1,"W":2,"ttbckg":3}

## training
hep_data_name = config.hep_data.split('/')[-1].split('.')[0]
output = "{0}/{1}".format( config.output_path,hep_data_name)
output += "/training-{0}/".format(date)
dnn.output_dir = output

vb.INFO("RUN :  Saving output to {0}".format(output))
if not os.path.isdir(output):
    vb.WARNING("RUN : '{0}' does not exist ".format(output))
    vb.WARNING("RUN :       Creating the directory. ")
    os.system( 'mkdir -p {0}'.format(output) )

## -- Copy the configuration file to the output directory
os.system("cp {0} {1}".format(sys.argv[1],output))

## -- Slice rows of the dataframe (remove from training)
##    list of strings with arguments separated by a space
slices = ['AK4_deepCSVb >= 0','AK4_deepCSVbb >= 0',
          'AK4_deepCSVc >= 0','AK4_deepCSVl >= 0']   # want all AK4 to have 'good' b-tagging scores

## Setup
dnn.initialize()

dnn.load_data(['target'])   # load HEP data (add 'target' branch to dataframe)
dnn.preprocess_data(slices) # equal statistics for each class & remove bad rows
dnn.train()                 # build and train the model!

## END ##
