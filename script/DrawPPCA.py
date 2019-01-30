import numpy as np
import matplotlib.pyplot as plt
import pylab
import seaborn as sns
import codecs

def main():
    draw_LL()
    draw_Z()

def draw_LL():
    file_name = "result/ppca/LLs.txt"
    lines = open(file_name, "r").readlines()
    LLs = list()
    for line in lines:
        LLs.append(float(line))
    fig = plt.figure()
    plt.plot(np.arange(len(LLs)), LLs)
    plt.xlabel("iteration")
    plt.ylabel("log-likelihood")
    plt.title("Transition of log-likelihood")
    fig.savefig("result/ppca/figure/LogLikelihood.png", dpi=300)
    plt.close(1)

def draw_Z():
    file_name = "result/ppca/Z.tsv"
    populations, colors, types = load_color_list()
    lines = open(file_name, "r").readlines()
    Z = np.zeros([len(lines), len(lines[0].split("\t"))-1])
    for i,line in enumerate(lines):
        temp_list = line.split("\t")
        for j,temp in enumerate(temp_list[:-1]):
            Z[i,j] = float(temp)
    n_components = len(Z[0])
    i_list = range(1,n_components,2)
    for i in i_list:
        left = list(); height = list()
        for j in range(len(Z)):
            left.append(Z[j,i-1])
            height.append(Z[j,i])
        fig = plt.figure()
        for type, color in zip(types, colors):
            temp_left = list()
            temp_height = list()
            for j,population in enumerate(populations):
                if(population == type):
                    temp_left.append(left[j])
                    temp_height.append(height[j])
            plt.scatter(temp_left, temp_height, c=color, label=type)
        plt.xlabel(str(i) + "-Component")
        plt.ylabel(str(i+1) + "-Component")
        plt.title("Probabilistic PCA")
        plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", borderaxespad=0)
        pylab.subplots_adjust(right=0.75)
        fig.savefig("result/ppca/figure/PPCA-" + str(int((i+1)/2)) + ".png",
                    dpi=300)
        plt.close(1)

def load_color_list():
    genotype_file = "data/genotype_reshaped.tsv"
    lines = open(genotype_file, "r").readlines()
    samples = lines[0].split("\t")[5:]
    population_file = "data/population.txt"
    lines = codecs.open(population_file, "r").readlines()
    IGRC_list = list(); NSFTV_list = list(); population_list = list()
    for line in lines:
        temp_list = line.split("\t")
        IGRC_list.append(temp_list[4])
        NSFTV_list.append(temp_list[5])
        population_list.append(temp_list[8])
    actual_populations = list()
    for sample in samples:
        for IGRC, NSFTV, population in zip(IGRC_list, NSFTV_list,
                population_list):
            if(IGRC + "@" in sample or NSFTV + "@" in sample):
                actual_populations.append(population)
    population_types = list()
    for population in actual_populations:
        if(population not in population_types):
            population_types.append(population)
    type_num = len(population_types)
    colors = sns.color_palette("hls", type_num)
    return actual_populations, colors, population_types

if __name__ == "__main__":
    main()
