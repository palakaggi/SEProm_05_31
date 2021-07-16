from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import *
from sklearn.model_selection import *
from mlxtend.feature_selection import SequentialFeatureSelector
from genetic_selection import GeneticSelectionCV
from sklearn.ensemble import RandomForestClassifier

def RecursiveFeatureSelection(df):
    x = df.drop('TSS', axis = 1)
    y = df['TSS']
    nof_list = np.arange(1,31)
    high_score = 0
    nof = 0
    score_list = []
    for n in range(len(nof_list)):
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=1)
        model1 = LogisticRegression()
        rfe1 = RFE(model1, nof_list[n])
        x_train_rfe1 = rfe1.fit_transform(x_train, y_train)
        x_test_rfe1 = rfe1.transform(x_test)
        model1.fit(x_train_rfe1, y_train)
        score = model1.score(x_test_rfe1, y_test)

        score_list.append(score)
        if score>high_score:
            high_score = score
            nof = nof_list[n]

        # model2 = RandomForestClassifier(n_estimators=50)
        # rfe2 = RFE(model2, nof_list[n])
        # x_train_rfe2 = rfe2.fit_transform(x_train, y_train)
        # x_test_rfe2 = rfe2.transform(x_test)
        # model2.fit(x_train_rfe2, y_train)
        # score = model2.score(x_test_rfe2, y_test)

        # score_list.append(score)
        # if score > high_score:
        #     high_score = score
        #     nof = nof_list[n]

    print("optimum number of features is", nof)
    print("score with %d features: %f" %(nof, high_score))

def forwardFeatureSelection(df):
    x = df.drop('TSS', axis = 1)
    y = df['TSS']
    x_train,x_test,y_train,y_test= train_test_split(x,y,test_size=0.2, random_state=1)

    # create the SequentialFeatureSelector object, and configure the parameters.
    sfs = SequentialFeatureSelector(RandomForestClassifier(),
                                    k_features=22,
                                    forward=True,
                                    floating=False,
                                    scoring='accuracy',
                                    cv=2)

    # fit the object to the training data.
    sfs = sfs.fit(x_train, y_train)

    # print the selected features.
    selected_features = x_train.columns[list(sfs.k_feature_idx_)]
    print(selected_features)

    # print the final prediction score.
    print(sfs.k_score_)

    # transform to the newly selected features.
    x_train_sfs = sfs.transform(x_train)
    x_test_sfs = sfs.transform(x_test)

def backwardFeatureSelection(df):
    # just set forward=False for backward feature selection.
    # create theSequentialFeatureSelector object, and configure the parameters.
    x = df.drop('TSS', axis=1)
    y = df['TSS']
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=1)
    sbs = SequentialFeatureSelector(LogisticRegression(),
                                    k_features=22,
                                    forward=False,
                                    floating=False,
                                    scoring='accuracy',
                                    cv=2)

    # fit the object to our training data.
    sbs = sbs.fit(x_train, y_train)

    # print the selected features.
    selected_features = x_train.columns[list(sbs.k_feature_idx_)]
    print(selected_features)

    # print the final prediction score.
    print(sbs.k_score_)

    # transform to the newly selected features.
    x_train_sfs = sbs.transform(x_train)
    x_test_sfs = sbs.transform(x_test)
