import pandas as pd
import numpy as np
import tensorflow as tf
import csv
import os
import glob
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, GaussianNoise, Reshape
from tensorflow.keras.layers import Flatten, LocallyConnected1D
from keras.backend.tensorflow_backend import set_session
from itertools import groupby


config = tf.ConfigProto()
# config.gpu_options.per_process_gpu_memory_fraction = 0.1
config.gpu_options.allow_growth = True
set_session(tf.Session(config=config))

# helper function
def split_text(s):
    for k, g in groupby(s, str.isalpha):
        yield ''.join(g)


# "build network"-function

def build_network(arc, drop_rate, dg, lc, x_train):
    def add_drops(model, drop_out, k):
        if dg[k].upper() == "D":
            model.add(Dropout(drop_out[0]))
        elif dg[k].upper() == "G":
            model.add(GaussianNoise(drop_out[k]))
        else:
            pass
        return model

    dg = dg.strip().split(",")
    arc = arc.strip().split(",")
    archit = []
    drop_rate = drop_rate.strip().split(",")
    drop_out_layers = []
    for layer in arc:
        archit.append(int(layer))
    for drops in drop_rate:
        drop_out_layers.append(float(drops))
    my_model = Sequential()
    if lc:
        my_model.add(Reshape(input_shape=(x_train.shape[1],), target_shape=(x_train.shape[1], 1)))
        my_model.add(LocallyConnected1D(1, 10, strides=7, input_shape=(x_train.shape[1], 1)))
        my_model.add(Flatten())
        start = 0
        my_model = add_drops(model=my_model, drop_out=drop_out_layers, k=start)
    else:
        my_model.add(Dense(archit[0], kernel_initializer=my_initializer, activation=my_act, input_shape=(x_train.
                                                                                                         shape[1],)))
        start = 1
        my_model = add_drops(model=my_model, drop_out=drop_out_layers, k=start)
    for i in range(start, len(archit)):
        my_model.add(Dense(archit[i], kernel_initializer=my_initializer, activation=my_act))
        my_model = add_drops(model=my_model, drop_out=drop_out_layers, k=i)
    my_model.add(Dense(expected_result_dimension, kernel_initializer=my_initializer))
    return my_model

def genopred(phenocsv, genocsv, cvfcsv):
    x = pd.read_csv(genocsv, engine="python")
    y = pd.read_csv(phenocsv, engine="python")
    cvfolds = pd.read_csv(cvfcsv, engine="python")
    y = y["phenotype_value"]
    rms = []
    for k in list(y.index):
        if np.isnan(y[k]):
            rms.append(k)
    for l in rms:
        y.pop(l)
        x.drop(l, inplace=True)
    for i in range(1, 6):
        print("Predicting cv_fold {} of {}".format(i, 5))
        y_train = []
        x_train = pd.DataFrame()
        y_test = []
        x_test = pd.DataFrame()
        for j in list(y.index):
            if cvfolds["cv_{}".format(i)][j] == 1:
                y_test.append(y[j])
                x_test = x_test.append(x.iloc[:, j+1])
            else:
                y_train.append(y[j])
                x_train = x_train.append(x.iloc[:,j+1])
        y_train = pd.DataFrame(y_train)
        x_train = pd.DataFrame(x_train)
        y_test = pd.DataFrame(y_test)
        x_test = pd.DataFrame(x_test)
        cv_model = build_network(arc=my_arc, drop_rate=my_drop_rate, dg=DG, lc=LC, x_train=x_train)
        cv_model.compile(loss=my_loss, optimizer=my_optimizer)
        cv_model.fit(x_train, y_train, epochs=my_epochs)
        y_hat = cv_model.predict(x_test)
        accuracy = np.corrcoef(y_hat.flatten(), np.asarray(y_test).flatten())
        foo = os.path.basename(phenocsv)
        mod = list(split_text(foo))[1]
        herit = list(split_text(foo))[0]
        result = [accuracy[0][1], mod, herit]
        print("\n===================================\nPrediction accuracy is at {}.".format(accuracy[0][1]))
        filename = "/home/s340454/master/genopredANN.csv"
        if os.path.isfile(filename):
            with open(filename, "a") as f:
                writer = csv.writer(f)
                writer.writerow(result)
        else:
            pd.DataFrame(accuracy).to_csv(filename)
        print("\nDone\n")


my_act = "sigmoid"
my_drop_rate = str('0.5')
my_arc = str('50,30')
DG = 'D,D,D,D,D,G'
LC = True
my_loss = "mean_squared_error"
my_optimizer = "Adam"
my_epochs = 100
my_initializer = 'truncated_normal'
expected_result_dimension = 1

genocsv = "/storage/full-share/gp_at/filtered_genos_csv/FT10.csv"
cvfcsv = "/storage/full-share/gp_at/cvf/FT10.csv"

# pipeline for final experiment (predicting all simluated phenotypes)
for phenocsv in glob.iglob("/home/s340454/master/sim_phenotypes/**", recursive=True):
    if os.path.isfile(phenocsv):
        genopred(phenocsv=phenocsv, genocsv=genocsv, cvfcsv=cvfcsv)
        
# pipeline for grid search of hyperparameters
#activations = ["relu", "sigmoid"]
#optimizers = ["Adam"]
#architectures =[str("50,30"), str("50,50"), str("50,30,15")]
#dropouts = [str("0.5,0.5,0.5")]
#phenocsv = "/storage/full-share/gp_at/new_phenos/FT10.csv"
#x = pd.read_csv(genocsv, engine="python")
#y = pd.read_csv(phenocsv, engine="python")
#cvfolds = pd.read_csv(cvfcsv, engine="python")
#y = y["phenotype_value"]
#rms = []
#for k in list(y.index):
#    if np.isnan(y[k]):
#        rms.append(k)
#    for l in rms:
#        y.pop(l)
#        x.drop(l, inplace=True)
#result = []
#
#for activation in activations:
#    for optimizer in optimizers:
#        for architecture in architectures:
#            for dropout in dropouts:
#                if dropout == "0.1,0.1,0.1":
#                    foo = 0.1
#                else:
#                    foo = 0.5
#                for i in range(1, 6):
#                    print("Predicting cv fold {} of {}".format(i, 5))
#                    y_train = []
#                    x_train = pd.DataFrame()
#                    y_test = []
#                    x_test = pd.DataFrame()
#                    for j in list(y.index):
#                        if cvfolds["cv_{}".format(i)][j] == 1:
#                            y_test.append(y[j])
#                            x_test = x_test.append(x.iloc[:, j+1])
#                        else:
#                            y_train.append(y[j])
#                            x_train = x_train.append(x.iloc[:, j+1])
#                    y_train = pd.DataFrame(y_train)
#                    y_test = pd.DataFrame(y_test)
#                    x_train = pd.DataFrame(x_train)
#                    x_test = pd.DataFrame(x_test)
#                    test_model = build_network(arc=architecture, drop_rate=dropout, dg=DG, lc=LC, x_train=x_train)
#                    test_model.compile(loss="MSE", optimizer=optimizer)
#                    test_model.fit(x_train, y_train, epochs=100)
#                    y_hat = test_model.predict(x_test)
#                    accuracy = np.corrcoef(y_hat.flatten(), np.asarray(y_test).flatten())
#                    print("\n==========================\nPrediction accuracy is at {}".format(accuracy[0][1]))
#                    result.append({"activation": activation, "optimizer": optimizer, "architecture": architecture, "dropout":foo, "accuracy":accuracy[0][1]})
#                    print("\nDone")
#
#final = pd.DataFrame.from_dict(result)