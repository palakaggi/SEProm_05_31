from readAndCreatingDF import readSequenceFile, getParameterDetails, createDataFrame
from plotting import plotting
import numpy as np
import matplotlib.pyplot as plt
import asyncPython
from multiprocessing import *

import datetime

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

    sequence_list_tss = list(sequence_map_tss.values())
#     sequence_list_cds = list(sequence_map_cds.values())

    tss_normalised_map = asyncPython.normalize_params(sequence_list_tss)

    print("normalising done")
    normalising_time = datetime.datetime.now()

#     cds_normalized_map = asyncPython.main(sequence_list_cds)
#     print(cds_normalized_map[0])

    struct_energy_map = asyncPython.energyStructParamsMP(tss_normalised_map)

    combining_params_time = datetime.datetime.now()

    createDataFrame.createMotifDF(struct_energy_map)
    createDataFrame.createDF(tss_normalised_map)

    df_time = datetime.datetime.now()

    print(reading_end - start)
    print(normalising_time - reading_end)
    print(combining_params_time-normalising_time)
    print(df_time-combining_params_time)


# #PLOTTING
# name = filepath_tss.split('/')[-1].split(' ')[0]
# print(name)
# # plotting(parameter_map_tss, parameter_map_cds, name)
# #PLOTTING OVER
