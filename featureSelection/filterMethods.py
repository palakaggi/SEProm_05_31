# import the required functions and object.
from sklearn.feature_selection import mutual_info_classif
from sklearn.feature_selection import SelectKBest
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn.feature_selection import f_classif

def mutualInformation(df):
    # select the number of features you want to retain.
    select_k = 15
    x = df.drop("TSS", axis=1)
    y = df["TSS"]
    x_train,x_test,y_train,y_test= train_test_split(x,y,test_size=0.2, random_state=1)

    # create the SelectKBest with the mutual info strategy.
    selection = SelectKBest(mutual_info_classif, k=select_k).fit(x_train, y_train)

    # display the retained features.
    features = x_train.columns[selection.get_support()]
    print(features)

def ANOVAFeatureSelection(df):
    select_k = 15
    x = df.drop("TSS", axis=1)
    y = df["TSS"]
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=1)

    # create the SelectKBest with the mutual info strategy.
    selection = SelectKBest(f_classif, k=select_k).fit(x_train, y_train)

    # display the retained features.
    features = x_train.columns[selection.get_support()]
    print(features)