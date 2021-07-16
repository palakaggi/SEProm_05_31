import numpy as np
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def corr_with_output1(combined_df):
    x = combined_df.drop("TSS", axis=1)
    y = combined_df["TSS"]
    corr_arr = []
    corr_dict = {}
# CALCULATING CORRELATION BETWEEN
    for i in x:
        corr1, _ = spearmanr(combined_df[i], y)
        corr_arr.append(corr1)
        corr_dict[i] = abs(corr1)
    # print(np.abs(corr_arr))
    return corr_dict
#CREATING CORRELATION BAR GRAPH
    # params = x.columns
    # plt.barh(params, corr_arr)
    # # plt.ylim((-1, 1))
    # print(len(corr_arr))
    # for index,value in enumerate(corr_arr):
    #     print(index, value)
    #     plt.text(value, index, str(value))
    # plt.show()


def correlation_pair(df):
    df_corr = df.drop("TSS", axis=1)
    cor_matrix = df_corr.corr(method = 'spearman').abs()
    # heatmap = sns.heatmap(cor_matrix, cmap = 'YlGnBu')
    # plt.show()
    # import sys
    # sys.exit()
    upper_tri = cor_matrix.where(np.triu(np.ones(cor_matrix.shape), k=1).astype(np.bool))
    sorted_mat = upper_tri.unstack().sort_values()
    reversed_sorted = sorted_mat.iloc[::-1].dropna()
    if reversed_sorted.iloc[0]>0.70:
        return reversed_sorted.index.values[0]
    else:
        return None
    # print(count)
    # import sys
    # sys.exit()
    # to_drop = [column for column in upper_tri.columns if any(upper_tri[column] > 0.6)]
    # print(to_drop)
    # df1 = df.drop(to_drop, axis=1)
    # print(df1.head())


def greedy_algo(combined_df, corr_with_op):
    # if len(combined_df.columns) == 15:
    #     return combined_df

    if 'motifs' in combined_df.columns:
        combined_df = combined_df.drop('motifs', axis=1)

    # corr_with_op = corr_with_output1(combined_df)

    # for i in range(len(combined_df.columns)-15):
    corr_pair = correlation_pair(combined_df)
    if corr_pair == None:
        return combined_df
    print(corr_pair)
    print(corr_with_op[corr_pair[0]],corr_with_op[corr_pair[1]])

    if corr_with_op[corr_pair[1]]<=corr_with_op[corr_pair[0]]:
        new_eliminated = combined_df.drop(corr_pair[1], axis = 1)
    else:
        new_eliminated = combined_df.drop(corr_pair[0], axis = 1)

    print(new_eliminated.columns)
    print(len(new_eliminated.columns))

    return greedy_algo(new_eliminated, corr_with_op)


