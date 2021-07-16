from readAndCreatingDF import readSequenceFile
import datetime
import pandas as pd
from multiprocessing import Pool


df = pd.read_excel("/Users/palakaggarwal/Desktop/Palak/Projects/SEProm_04_28/train_data/triAndTetraData/Trimer.xlsx",sheet_name="Sheet3",index_col=0)

start = datetime.datetime.now()
trimer_tss = '/Users/palakaggarwal/Desktop/Palak/Projects/SEProm_04_28/train_data/tss.txt'
trimer_cds = '/Users/palakaggarwal/Desktop/Palak/Projects/SEProm_04_28/train_data/cds.txt'

# sequence_map_tss = readSequenceFile.readsequencefile(path_tss)
# print(sequence_map_tss)

# sequence_map_cds = readSequenceFile.readsequencefile(path_cds)
# print(sequence_map_cds)


def iterateSequences(seq_map,dataframe):
    """
     :param sequence_map:
     :return: returns a list of dicts for each sequence
     """
    sequence_list = list(seq_map.values())
    print("starting normalisation of read sequences")
    pool = Pool()
    param_list = pool.starmap(calculateParameters, [(seq, df) for seq in sequence_list])
    pool.close()
    return param_list


def calculateParameters(sequence,dframe):
    param_map = {"a":[], "b":[], "c":[], "d":[], "f":[], "g":[], "h":[], "i": [], "j": [], "k":[], "t":[], "u":[], "v":[], "w":[], "x":[], "y":[], "z":[], "aa":[], "ab":[], "ac":[], "ad":[], "ae":[] }
    noofbases = len(sequence)

    if noofbases == 0:
        return
    trimotifs=[]
    for m in range(noofbases-2):
        trimotifs.append(sequence[m:m+3])
    for motif in trimotifs:
        # m=str(motif, 'utf-8')
        if (motif in dframe.columns):
            assign_params(param_map,dframe[motif])

    return calculateMovingAverages(param_map)

def assign_params(param_map, plist):
    param_map['a'].append(plist['a'])
    param_map['b'].append(plist["b"])
    param_map['c'].append(plist["c"])
    param_map['d'].append(plist["d"])
    param_map['f'].append(plist["f"])
    param_map['g'].append(plist["g"])
    param_map['h'].append(plist["h"])
    param_map['i'].append(plist["i"])
    param_map['j'].append(plist["j"])
    param_map['k'].append(plist["k"])
    param_map['t'].append(plist["t"])
    param_map['u'].append(plist["u"])
    param_map['v'].append(plist["v"])
    param_map['w'].append(plist["w"])
    param_map['x'].append(plist["x"])
    param_map['y'].append(plist["y"])
    param_map['z'].append(plist["z"])
    param_map['aa'].append(plist["aa"])
    param_map['ab'].append(plist["ab"])
    param_map['ac'].append(plist["ac"])
    param_map['ad'].append(plist["ad"])
    param_map['ae'].append(plist["ae"])
    return param_map

def calculateMovingAverages(param_map):
    moving_win_size = 25
    moving_param_map = {}
    for k, v in param_map.items():
        arr = v
        moving_param_map[k] = []
        for i in range(0, len(arr) - moving_win_size + 1):
            sum = 0
            for j in range(i, i + moving_win_size):
                sum += arr[j]
            avg = sum / moving_win_size
            moving_param_map[k].append(avg)

    return moving_param_map


if __name__ == '__main__':
    sequence_map_tss = readSequenceFile.readSequenceFile(trimer_tss)
    sequence_map_cds = readSequenceFile.readSequenceFile(trimer_cds)

    print("Reading of trimer data done!!")

    tri_tss_seq_list = iterateSequences(sequence_map_tss, df)
    tetra_cds_seq_list = iterateSequences(sequence_map_cds, df)

    print("Normalization of trimer data done!!")

    #plotting
    """
    name = tetramer_tss.split('/')[-1].split(' ')[0]
    plotting.plotting_tetra(tetra_tss_seq_list, tetra_cds_seq_list, name)
    print("Plotting of Tetramer data done!!")
    """
