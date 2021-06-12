import pandas as pd
import multiprocessing as mp

params = ['a','b','c','d','e','f','g','h','i','j','k','l','ma','n','o','p','q','r','s','t','u','v','w','x','y','z','aa','ab','ac','ad','ae']
motifParams = ['structuralIncreasing_params','structuralDecreasing_params','energyIncreasing_params','energyDecreasing_params']


def avg_window(start, stop, params, params_map):
    arr = []
    for p in params:
        try:
            params_list = list(params_map[p].values())
        except:
            params_list = params_map[p]
        avg_p = sum(params_list[start:stop])/(stop-start)
        arr.append(float(avg_p))
    return arr

def createMotifDF(normalized_map_tss):
    """
    :param normalized_map_tss: Map containing 976 values per seq per paramater
    :return: A pandas DF with columns(SI, SD, EI, ED and TSS, Motif)
    """
    print("CREATING MOTIF DATAFRAME")
    seq_data = pd.DataFrame(columns=motifParams)
    # print(seq_data)
    params1 = motifParams+['motif']
    # ADDING TSS DATA TO DATAFRAME
    for seq in normalized_map_tss.keys():
        y = len(seq_data.index)
        m= avg_window(486, 493, motifParams, normalized_map_tss[seq])
        m.append('m_0_0')
        # m = pd.Series(m,index = params1)
        n = avg_window(462, 472, motifParams, normalized_map_tss[seq])
        n.append('m_0_1')
        l = avg_window(420, 437, motifParams, normalized_map_tss[seq])
        l.append('m_0_2')
        seq_data = seq_data.append(pd.Series(m,index = params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(n,index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
        # print(seq_data)
        # import sys
        # sys.exit()
        m = avg_window(487, 494, motifParams, normalized_map_tss[seq])
        m.append('m_1_0')

        n = avg_window(464, 471, motifParams, normalized_map_tss[seq])
        n.append('m_1_1')
        l = avg_window(453, 462, motifParams, normalized_map_tss[seq])
        l.append('m_1_2')
        k = avg_window(440, 452, motifParams, normalized_map_tss[seq])
        k.append('m_1_3')
        seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

        m = avg_window(483, 494, motifParams, normalized_map_tss[seq])
        m.append('m_2_0')
        n = avg_window(463, 472, motifParams, normalized_map_tss[seq])
        n.append('m_2_1')
        l = avg_window(450, 459, motifParams, normalized_map_tss[seq])
        l.append('m_2_2')
        k = avg_window(434, 445, motifParams, normalized_map_tss[seq])
        k.append('m_2_3')
        seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

        m = avg_window(487, 498, motifParams, normalized_map_tss[seq])
        m.append('m_3_0')
        n = avg_window(463, 472, motifParams, normalized_map_tss[seq])
        n.append('m_3_1')
        l = avg_window(450, 459, motifParams, normalized_map_tss[seq])
        l.append('m_3_2')
        k = avg_window(432, 445, motifParams, normalized_map_tss[seq])
        k.append('m_3_3')
        seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

        z = len(seq_data.index)

        m = avg_window(786, 793, motifParams, normalized_map_tss[seq])
        m.append('nm_0_0')
        n = avg_window(763, 772, motifParams, normalized_map_tss[seq])
        n.append('nm_0_1')
        l = avg_window(720, 737, motifParams, normalized_map_tss[seq])
        l.append('nm_0_2')
        seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)

        m = avg_window(787, 794, motifParams, normalized_map_tss[seq])
        m.append('nm_1_0')
        n = avg_window(764, 771, motifParams, normalized_map_tss[seq])
        n.append('nm_1_1')
        l = avg_window(753, 762, motifParams, normalized_map_tss[seq])
        l.append('nm_1_2')
        k = avg_window(740, 752, motifParams, normalized_map_tss[seq])
        k.append('nm_1_3')
        seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

        m = avg_window(783, 794, motifParams, normalized_map_tss[seq])
        m.append('nm_2_0')
        n = avg_window(763, 772, motifParams, normalized_map_tss[seq])
        n.append('nm_2_1')
        l = avg_window(750, 759, motifParams, normalized_map_tss[seq])
        l.append('nm_2_2')
        k = avg_window(734, 745, motifParams, normalized_map_tss[seq])
        k.append('nm_2_3')
        seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

        m = avg_window(787, 798, motifParams, normalized_map_tss[seq])
        m.append('nm_3_0')
        n = avg_window(763, 772, motifParams, normalized_map_tss[seq])
        n.append('nm_3_1')
        l = avg_window(750, 759, motifParams, normalized_map_tss[seq])
        l.append('nm_3_2')
        k = avg_window(732, 745, motifParams, normalized_map_tss[seq])
        k.append('nm_3_3')
        seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
        seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

        seq_data.loc[y:z, 'TSS'] = int(1)
        seq_data.loc[z:, 'TSS'] = int(0)



    # seq_data['SI'] = seq_data['k']+seq_data['p']+seq_data['q']+seq_data['s']+seq_data['u']+seq_data['v']+seq_data['x']+seq_data['g']+seq_data['n']+seq_data['r']+seq_data['aa']+seq_data['ab']
    # seq_data = seq_data.drop(['k','p','q','s','u','v','x','g','n','r','aa','ab'], axis=1)
    # seq_data['SD'] = seq_data['a']+seq_data['b']+seq_data['f']+seq_data['h']+seq_data['l']+seq_data['ma']+seq_data['c']+seq_data['d']+seq_data['e']+seq_data['i']+seq_data['j']+seq_data['o']+seq_data['t']+seq_data['w']+seq_data['y']+seq_data['z']
    # seq_data = seq_data.drop(['a','b','f','h','l','ma','c','d','e','i','j','o','t','w','y','z'], axis = 1)
    # seq_data['EI'] = seq_data['ac']+seq_data['ad']
    # seq_data = seq_data.drop(['ac','ad'], axis=1)
    # seq_data['ED'] = seq_data['ae']
    # seq_data = seq_data.drop('ae', axis=1)

    seq_data.to_csv('TEST_MOTIF_DATA.csv')
    return

def createDF(normalized_map_tss):
    """
    :param normalized_map_tss: Map containing 976 values per seq per paramater{0:a:[],b:[],c:[]...,
                                                                                1:a:[],b:[]...}
    :return: A pandas DF with 31 parameters as their columns and TSS/no TSS
    """
    print("CREATING NORMALISED DATAFRAME")
    nml = len(normalized_map_tss)
    pool = mp.Pool(4)
    m_list = pool.starmap(avg_window,[(425, 505, params, normalized_map_tss[seq]) for seq in range(nml)])
    pool.close()

    # print(m_list)
    seq_data = pd.DataFrame(m_list, columns=params)
    seq_data['TSS'] = 1
    # print(seq_data)
    # for seq in normalized_map_tss.keys():
    #     m = [float(i) for i in avg_window(425, 505, params, normalized_map_tss[seq])]
    #     seq_data.loc[seq] = m

    # l = len(seq_data.index)

    pool = mp.Pool(4)
    m_list = pool.starmap(avg_window, [(700, 780, params, normalized_map_tss[seq]) for seq in range(nml)])
    pool.close()

    # print(m_list)
    seq_data_no_tss = pd.DataFrame(m_list, columns=params)
    seq_data_no_tss['TSS'] = 0

    # print(seq_data_no_tss)

    final_df = pd.concat([seq_data, seq_data_no_tss], ignore_index=True)
    # print(final_df)

    return
    # for seq in normalized_map_tss.keys():
    #     m = [float(i) for i in avg_window(700, 780, params, normalized_map_tss[seq])]
    #     seq_data.loc[l + seq] = m

    # seq_data.loc[:l, 'TSS'] = 1
    # seq_data.loc[l:, 'TSS'] = 0

    # tss['TSS'] = 1
    # no_tss['TSS'] = 0

    # combined_df = pd.concat([tss, no_tss], ignore_index=True)
    # seq_data.to_csv('training80window_no_mov_avg.csv')



