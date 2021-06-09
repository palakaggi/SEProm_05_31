from readAndCreatingDF import readSequenceFile, getParameterDetails, createDataFrame
from plotting import plotting
import numpy as np
import matplotlib.pyplot as plt
import asyncPython
from multiprocessing import *

import datetime

# print(datetime.datetime.now())

# def main(sequence_map):
if __name__ == '__main__':
    start=datetime.datetime.now()

    filepath_tss = "/Users/palakaggarwal/Desktop/Palak/Projects/SEProm_04_28/train_data/tss.txt"
    filepath_cds = "/Users/palakaggarwal/Desktop/Palak/Projects/SEProm_04_28/train_data/cds.txt"

    try:
        f = open(filepath_tss)
    except NameError:
        filepath_tss = str(input("Please enter the input sequence file."))

    sequence_map_tss = readSequenceFile.readSequenceFile(filepath_tss)
    sequence_map_cds = readSequenceFile.readSequenceFile(filepath_cds)
    print("reading done")
    reading_end = datetime.datetime.now()


    sequence_list = list(sequence_map_tss.values())
#     print(type(sequence_list[0]))

    asyncPython.main(sequence_list)
    print(reading_end - start)
    print(datetime.datetime.now()-reading_end)

    import sys
    sys.exit()
# print(sequence_map_tss.keys())


# parameter_map_tss = {}
# parameter_map_tss = getParameterDetails.iterateSequences(sequence_map_tss)
#
# parameter_map_cds = {}
# parameter_map_cds = getParameterDetails.iterateSequences(sequence_map_cds)
#
#
# # print(parameter_map_tss['combined_params_map'][0])
#
# #PLOTTING
# name = filepath_tss.split('/')[-1].split(' ')[0]
# print(name)
# # plotting(parameter_map_tss, parameter_map_cds, name)
# #PLOTTING OVER
#
# #CREATE DATAFRAME
# createDataFrame.createMotifDF(parameter_map_tss['combined_params_map'])
# createDataFrame.createDF(parameter_map_tss['normalized_params_map'])
#
# print(datetime.datetime.now())
#
