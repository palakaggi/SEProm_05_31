import pandas as pd


params = ['a','b','c','d','e','f','g','h','i','j','k','l','ma','n','o','p','q','r','s','t','u','v','w','x','y','z','aa','ab','ac','ad','ae']

def avg_window(start, stop, params, params_map):
    arr = []
    for p in params:
        avg_p = sum(params_map[p][start:stop])/(stop-start)
        arr.append(avg_p)
    return arr

def createMotifDF(normalized_map_tss):
    """
    :param normalized_map_tss: Map containing 976 values per seq per paramater
    :return: A pandas DF with columns(SI, SD, EI, ED and TSS, Motif)
    """
    print("CREATING MOTIF DATAFRAME")
    seq_data = pd.DataFrame(columns=params)
    # print(seq_data)
    params1 = ['a','b','c','d','e','f','g','h','i','j','k','l','ma','n','o','p','q','r','s','t','u','v','w','x','y','z','aa','ab','ac','ad','ae', 'motifs']
    # ADDING TSS DATA TO DATAFRAME
    for seq in normalized_map_tss.keys():
        y = len(seq_data.index)
        m = pd.Series(avg_window(486,493,params,normalized_map_tss[seq]).append('m_0_0'),index = params1)
        n = pd.Series(avg_window(462, 472, params, normalized_map_tss[seq]).append('m_0_1'), index=params1)
        l = pd.Series(avg_window(420, 437, params, normalized_map_tss[seq]).append('m_0_2'), index=params1)
        seq_data = seq_data.append(m, ignore_index=True)
        seq_data = seq_data.append(n, ignore_index=True)
        seq_data = seq_data.append(l, ignore_index=True)
        m = pd.Series(avg_window(487, 494, params, normalized_map_tss[seq]).append('m_1_0'), index=params1)
        n = pd.Series(avg_window(464, 471, params, normalized_map_tss[seq]).append('m_1_1'), index=params1)
        l = pd.Series(avg_window(453, 462, params, normalized_map_tss[seq]).append('m_1_2'), index=params1)
        k = pd.Series(avg_window(440, 452, params, normalized_map_tss[seq]).append('m_1_3'), index=params1)
        seq_data = seq_data.append(m, ignore_index=True)
        seq_data = seq_data.append(n, ignore_index=True)
        seq_data = seq_data.append(l, ignore_index=True)
        seq_data = seq_data.append(k, ignore_index=True)
        m = pd.Series(avg_window(483, 494, params, normalized_map_tss[seq]).append('m_2_0'), index=params1)
        n = pd.Series(avg_window(463, 472, params, normalized_map_tss[seq]).append('m_2_1'), index=params1)
        l = pd.Series(avg_window(450, 459, params, normalized_map_tss[seq]).append('m_2_2'), index=params1)
        k = pd.Series(avg_window(434, 445, params, normalized_map_tss[seq]).append('m_2_3'), index=params1)
        seq_data = seq_data.append(m, ignore_index=True)
        seq_data = seq_data.append(n, ignore_index=True)
        seq_data = seq_data.append(l, ignore_index=True)
        seq_data = seq_data.append(k, ignore_index=True)
        m = pd.Series(avg_window(487, 498, params, normalized_map_tss[seq]).append('m_3_0'), index=params1)
        n = pd.Series(avg_window(463, 472, params, normalized_map_tss[seq]).append('m_3_1'), index=params1)
        l = pd.Series(avg_window(450, 459, params, normalized_map_tss[seq]).append('m_3_2'), index=params1)
        k = pd.Series(avg_window(432, 445, params, normalized_map_tss[seq]).append('m_3_3'), index=params1)
        seq_data = seq_data.append(m, ignore_index=True)
        seq_data = seq_data.append(n, ignore_index=True)
        seq_data = seq_data.append(l, ignore_index=True)
        seq_data = seq_data.append(k, ignore_index=True)

        z = len(seq_data.index)

        m = pd.Series(avg_window(786, 793, params, normalized_map_tss[seq]).append('nm_0_0'), index=params1)
        n = pd.Series(avg_window(763, 772, params, normalized_map_tss[seq]).append('nm_0_1'), index=params1)
        l = pd.Series(avg_window(720, 737, params, normalized_map_tss[seq]).append('nm_0_2'), index=params1)
        seq_data = seq_data.append(m, ignore_index=True)
        seq_data = seq_data.append(n, ignore_index=True)
        seq_data = seq_data.append(l, ignore_index=True)
        m = pd.Series(avg_window(787, 794, params, normalized_map_tss[seq]).append('nm_1_0'), index=params1)
        n = pd.Series(avg_window(764, 771, params, normalized_map_tss[seq]).append('nm_1_1'), index=params1)
        l = pd.Series(avg_window(753, 762, params, normalized_map_tss[seq]).append('nm_1_2'), index=params1)
        k = pd.Series(avg_window(740, 752, params, normalized_map_tss[seq]).append('nm_1_3'), index=params1)
        seq_data = seq_data.append(m, ignore_index=True)
        seq_data = seq_data.append(n, ignore_index=True)
        seq_data = seq_data.append(l, ignore_index=True)
        seq_data = seq_data.append(k, ignore_index=True)
        m = pd.Series(avg_window(783, 794, params, normalized_map_tss[seq]).append('nm_2_0'), index=params1)
        n = pd.Series(avg_window(763, 772, params, normalized_map_tss[seq]).append('nm_2_1'), index=params1)
        l = pd.Series(avg_window(750, 759, params, normalized_map_tss[seq]).append('nm_2_2'), index=params1)
        k = pd.Series(avg_window(734, 745, params, normalized_map_tss[seq]).append('nm_2_3'), index=params1)
        seq_data = seq_data.append(m, ignore_index=True)
        seq_data = seq_data.append(n, ignore_index=True)
        seq_data = seq_data.append(l, ignore_index=True)
        seq_data = seq_data.append(k, ignore_index=True)
        m = pd.Series(avg_window(787, 798, params, normalized_map_tss[seq]).append('nm_3_0'), index=params1)
        n = pd.Series(avg_window(763, 772, params, normalized_map_tss[seq]).append('nm_3_1'), index=params1)
        l = pd.Series(avg_window(750, 759, params, normalized_map_tss[seq]).append('nm_3_2'), index=params1)
        k = pd.Series(avg_window(732, 745, params, normalized_map_tss[seq]).append('nm_3_3'), index=params1)
        seq_data = seq_data.append(m, ignore_index=True)
        seq_data = seq_data.append(n, ignore_index=True)
        seq_data = seq_data.append(l, ignore_index=True)
        seq_data = seq_data.append(k, ignore_index=True)

        seq_data.loc[y:z, 'TSS'] = int(1)
        seq_data.loc[z:, 'TSS'] = int(0)

    seq_data['SI'] = seq_data['k']+seq_data['p']+seq_data['q']+seq_data['s']+seq_data['u']+seq_data['v']+seq_data['x']+seq_data['g']+seq_data['n']+seq_data['r']+seq_data['aa']+seq_data['ab']
    seq_data = seq_data.drop(['k','p','q','s','u','v','x','g','n','r','aa','ab'], axis=1)
    seq_data['SD'] = seq_data['a']+seq_data['b']+seq_data['f']+seq_data['h']+seq_data['l']+seq_data['ma']+seq_data['c']+seq_data['d']+seq_data['e']+seq_data['i']+seq_data['j']+seq_data['o']+seq_data['t']+seq_data['w']+seq_data['y']+seq_data['z']
    seq_data = seq_data.drop(['a','b','f','h','l','ma','c','d','e','i','j','o','t','w','y','z'], axis = 1)
    seq_data['EI'] = seq_data['ac']+seq_data['ad']
    seq_data = seq_data.drop(['ac','ad'], axis=1)
    seq_data['ED'] = seq_data['ae']
    seq_data = seq_data.drop('ae', axis=1)

    seq_data.to_csv('TEST_MOTIF_DATA.csv')
    return

def createDF(normalized_map_tss):
    """
    :param normalized_map_tss: Map containing 976 values per seq per paramater{0:a:[],b:[],c:[]...,
                                                                                1:a:[],b:[]...}
    :return: A pandas DF with 31 parameters as their columns and TSS/no TSS
    """
    print("CREATING NORMALISED DATAFRAME")
    seq_data = pd.DataFrame(columns=params)
    for seq in normalized_map_tss.keys():
        m = [float(i) for i in avg_window(425, 505, params, normalized_map_tss[seq])]
        seq_data.loc[seq] = m

    l = len(seq_data.index)

    for seq in normalized_map_tss.keys():
        m = [float(i) for i in avg_window(700, 780, params, normalized_map_tss[seq])]
        seq_data.loc[l + seq] = m

    seq_data.loc[:l, 'TSS'] = 1
    seq_data.loc[l:, 'TSS'] = 0

    # tss['TSS'] = 1
    # no_tss['TSS'] = 0

    # combined_df = pd.concat([tss, no_tss], ignore_index=True)
    seq_data.to_csv('training80window_no_mov_avg.csv')



