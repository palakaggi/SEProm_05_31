from readAndCreatingDF import readSequenceFile
import pandas as pd
from plotting import plotting
from multiprocessing import Pool

tetramer_tss = '/Users/palakaggarwal/Desktop/Palak/Projects/SEProm_04_28/train_data/tss.txt'
paramVals = '/Users/palakaggarwal/Desktop/Palak/Projects/SEProm_04_28/train_data/triAndTetraData/Tetramer.xlsx'
tetramer_cds = '/Users/palakaggarwal/Desktop/Palak/Projects/SEProm_04_28/train_data/cds.txt'

strInc = []
strDec = []

def dataCleaning(df):
    # df = df.drop(['Tetramer'], axis =1)
    duplicates = df.duplicated(subset = ['l','m','n','o','p','q'], keep = False)
    df2 = df[~duplicates]
    df1 = df[duplicates]
    df1 = (df1.groupby(df1.columns.tolist())
           .apply(lambda x: tuple(x.index))
           .reset_index(name='Tetramer'))
    df1=df1.set_index('Tetramer')
    df = pd.concat([df1,df2])
    df = df.T
    return df

paramValues = pd.read_excel(paramVals, sheet_name='Sheet3', index_col=0)


def calculateParameters(sequence_map_per_seq,paramValues):
    param_map = {'l':[],'m':[], 'n':[], 'o':[], 'p':[], 'q':[]}
    # shift = slide= rise = tilt = roll = twist =0
    no_of_bases = len(sequence_map_per_seq)
    list_motifs = []
    for m in range(no_of_bases-3):
        list_motifs.append(sequence_map_per_seq[m:m+4])
    for motif in list_motifs:
        if motif in paramValues.columns:
            param_map['l'].append(paramValues[motif]['l'])
            param_map['m'].append(paramValues[motif]['m'])
            param_map['n'].append(paramValues[motif]['n'])
            param_map['o'].append(paramValues[motif]['o'])
            param_map['p'].append(paramValues[motif]['p'])
            param_map['q'].append(paramValues[motif]['q'])
        else:
            print("not found", motif)

    return calculateMovingAverages(param_map)

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

def iterateSequences(sequence_map):
    """
    :param sequence_map:
    :return: returns a list of dicts for each sequence
    """
    sequence_list = list(sequence_map.values())
    print("starting normalisation of read sequences")
    pool = Pool()
    param_list = pool.starmap(calculateParameters, [(seq, paramValues) for seq in sequence_list])
    pool.close()
    return param_list


if __name__ == '__main__':
    sequence_map_tss = readSequenceFile.readSequenceFile(tetramer_tss)
    # sequence_map_cds =readSequenceFile.readSequenceFile(tetramer_cds)

    print("Reading of tetramer data done!!")

    tetra_tss_seq_list = iterateSequences(sequence_map_tss)
    # tetra_cds_seq_list = iterateSequences(sequence_map_cds)

    print("Normalization of tetramer data done!!")

    #plotting
    """
    name = tetramer_tss.split('/')[-1].split(' ')[0]
    plotting.plotting_tetra(tetra_tss_seq_list, tetra_cds_seq_list, name)
    print("Plotting of Tetramer data done!!")
    """