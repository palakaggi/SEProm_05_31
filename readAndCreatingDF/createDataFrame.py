import pandas as pd
import multiprocessing as mp

params = ['a','b','c','d','e','f','g','h','i','j','k','l','ma','n','o','p','q','r','s','t','u','v','w','x','y','z','aa','ab','ac','ad','ae']
motifParams = ['SIParams_all_seq','SDParams_all_seq', 'EIparams_all_seq', 'EDParams_all_seq']

def motif_helper(normalised_seq_list):
    seq_data = pd.DataFrame(columns=motifParams)

    params1 = motifParams+['motif']
    m = avg_window(486, 493, motifParams, normalised_seq_list)
    m.append('m_0_0')
    n = avg_window(462, 472, motifParams, normalised_seq_list)
    n.append('m_0_1')
    l = avg_window(420, 437, motifParams, normalised_seq_list)
    l.append('m_0_2')
    seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
    # print(seq_data)

    m = avg_window(487, 494, motifParams, normalised_seq_list)
    m.append('m_1_0')
    n = avg_window(464, 471, motifParams, normalised_seq_list)
    n.append('m_1_1')
    l = avg_window(453, 462, motifParams, normalised_seq_list)
    l.append('m_1_2')
    k = avg_window(440, 452, motifParams, normalised_seq_list)
    k.append('m_1_3')
    seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

    m = avg_window(483, 494, motifParams, normalised_seq_list)
    m.append('m_2_0')
    n = avg_window(463, 472, motifParams, normalised_seq_list)
    n.append('m_2_1')
    l = avg_window(450, 459, motifParams, normalised_seq_list)
    l.append('m_2_2')
    k = avg_window(434, 445, motifParams, normalised_seq_list)
    k.append('m_2_3')
    seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

    m = avg_window(487, 498, motifParams, normalised_seq_list)
    m.append('m_3_0')
    n = avg_window(463, 472, motifParams, normalised_seq_list)
    n.append('m_3_1')
    l = avg_window(450, 459, motifParams, normalised_seq_list)
    l.append('m_3_2')
    k = avg_window(432, 445, motifParams, normalised_seq_list)
    k.append('m_3_3')
    seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

    m = avg_window(786, 793, motifParams, normalised_seq_list)
    m.append('nm_0_0')
    n = avg_window(763, 772, motifParams, normalised_seq_list)
    n.append('nm_0_1')
    l = avg_window(720, 737, motifParams, normalised_seq_list)
    l.append('nm_0_2')
    seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)

    m = avg_window(787, 794, motifParams, normalised_seq_list)
    m.append('nm_1_0')
    n = avg_window(764, 771, motifParams, normalised_seq_list)
    n.append('nm_1_1')
    l = avg_window(753, 762, motifParams, normalised_seq_list)
    l.append('nm_1_2')
    k = avg_window(740, 752, motifParams, normalised_seq_list)
    k.append('nm_1_3')
    seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

    m = avg_window(783, 794, motifParams, normalised_seq_list)
    m.append('nm_2_0')
    n = avg_window(763, 772, motifParams, normalised_seq_list)
    n.append('nm_2_1')
    l = avg_window(750, 759, motifParams, normalised_seq_list)
    l.append('nm_2_2')
    k = avg_window(734, 745, motifParams, normalised_seq_list)
    k.append('nm_2_3')
    seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)

    m = avg_window(787, 798, motifParams, normalised_seq_list)
    m.append('nm_3_0')
    n = avg_window(763, 772, motifParams, normalised_seq_list)
    n.append('nm_3_1')
    l = avg_window(750, 759, motifParams, normalised_seq_list)
    l.append('nm_3_2')
    k = avg_window(732, 745, motifParams, normalised_seq_list)
    k.append('nm_3_3')
    seq_data = seq_data.append(pd.Series(m, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(n, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(l, index=params1), ignore_index=True)
    seq_data = seq_data.append(pd.Series(k, index=params1), ignore_index=True)
    return seq_data

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
    print("-------CREATING MOTIF DATAFRAME-------")
    combined_params_list = list(normalized_map_tss.values())
    # print(len(combined_params_list))
    # print(combined_params_list[0])
    # print(type(combined_params_list[0]))
    pool = mp.Pool()
    results = pool.starmap(motif_helper, [[combined_params_list[seq]] for seq in range(len(combined_params_list))])
    pool.close()
    pool.join()
    seq_data = pd.concat(results)

    print(seq_data)

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

    seq_data_no_tss = pd.DataFrame(m_list, columns=params)
    seq_data_no_tss['TSS'] = 0

    final_df = pd.concat([seq_data, seq_data_no_tss], ignore_index=True)
    final_df.to_csv('training80window_no_mov_avg.csv')
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
