import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt
import seaborn as sns
import pylab

def linear_regression(g_training, g_test, p_training, p_test, pop_test, CV):
    num_ana = str(CV+1)
    X = g_training
    Y = p_training
    for exp_type in ["neutral", "lasso", "ridge"]:
        if(exp_type == "neutral"):
            clf = linear_model.LinearRegression()
        elif(exp_type == "lasso"):
            clf = linear_model.Lasso()
        elif(exp_type == "ridge"):
            clf = linear_model.Ridge()
        clf.fit(X,Y)
        weights = clf.coef_
        output = open("result/LinearRegression/weight_" + exp_type + "_" +\
                      num_ana + ".txt", "w")
        for weight in weights:
            output.write(str(weight) + "\n")
        prediction(g_test, p_test, pop_test, clf, CV, exp_type)
    

def prediction(g_test, p_test, pop_test, clf, CV, exp_type):
    X = g_test
    p_predicted = clf.predict(X)
    corrcoef = np.corrcoef(p_test, p_predicted)[0,1]
    pop_types = list()
    for pop in pop_test:
        if(pop not in pop_types):
            pop_types.append(pop)
    type_num = len(pop_types); colors = sns.color_palette("hls", type_num)
    fig = plt.figure()
    for pop_type, color in zip(pop_types, colors):
        temp_left = list()
        temp_height = list()
        for j,pop in enumerate(pop_test):
            if(pop == pop_type):
                temp_left.append(p_test[j])
                temp_height.append(p_predicted[j])
        plt.scatter(temp_left, temp_height, c=color, label=pop_type)
    plt.xlabel("Correct value")
    plt.ylabel("Predicted value")
    plt.title("Grain length (Correlation coefficient: " + str(corrcoef)+ ")")
    plt.legend(bbox_to_anchor=(1.01,1), loc="upper left", borderaxespad=0)
    pylab.subplots_adjust(right=0.75)
    name = "result/LinearRegression/figure/plot_" + exp_type + "_" +\
           str(CV+1) + ".png"
    plt.savefig(name, dpi=300)
    plt.close(1)
