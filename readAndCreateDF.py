from readAndCreatingDF import readSequenceFile, getParameterDetails, createDataFrame
from plotting import plotting
import numpy as np
import matplotlib.pyplot as plt

import datetime

print(datetime.datetime.now())

filepath_tss = "/Users/palakaggarwal/Desktop/dinucleotide_project/train_data/TSS_Seq/Salmonella enteric serovar typhimurium"
filepath_cds = "//Users/palakaggarwal/Desktop/dinucleotide_project/train_data/CDS_Seq/Salmonella_enteric_serovar_typhimurium_CDS"

try:
    f = open(filepath_tss)
except NameError:
    filepath_tss = str(input("Please enter the input sequence file."))

sequence_map_tss = readSequenceFile.readSequenceFile(filepath_tss)
sequence_map_cds = readSequenceFile.readSequenceFile(filepath_cds)

# print(sequence_map_tss.keys())

parameter_map_tss = {}
parameter_map_tss = getParameterDetails.iterateSequences(sequence_map_tss)

parameter_map_cds = {}
parameter_map_cds = getParameterDetails.iterateSequences(sequence_map_cds)

#PLOTTING
name = filepath_tss.split('/')
name = name[-1].split(' ')[0]
# plotting(parameter_map_tss, parameter_map_cds, name)
#PLOTTING OVER

#CREATE DATAFRAME
createDataFrame.createMotifDF(parameter_map_tss['normalized_params_map'])

p1rint(datetime.datetime.now())

