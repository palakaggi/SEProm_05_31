from statsmodels.stats.outliers_influence import variance_inflation_factor
import pandas as pd
import numpy as np

def cal_vif(df):
    x = df.drop('TSS', axis=1)
    y = df['TSS']

    thresh = 100
    output = pd.DataFrame()
    k = x.shape[1]

    vif = [variance_inflation_factor(x.values,i) for i in range(x.shape[1])]
    for i in range(1,k):
        print('Iteration number:', i)
        print(vif)
        a = np.argmax(vif)
        print('Max vif is for variable number', a+1)
        if vif[a]<=thresh:
            break
        if(i>=1):
            x = x.drop(x.columns[a], axis = 1)
            vif = [variance_inflation_factor(x.values,j) for j in range(x.shape[1])]
        # elif(i>1):
        #     output = output.drop(output.columns[a], axis=1)
        #     vif = [variance_inflation_factor(output.values, j) for j in range(output.shape[1])]
    return x.columns
    # print(k)
    # print(vif)