import numpy as np
import random
import LinearRegression

def main():
    markers, samples, genotype = get_genotype()
    phenotype_lib = get_phenotype(samples)
    population_lib = get_population(samples)
    num_samples = len(genotype); num_markers = len(genotype[0])
    sampling_list = range(num_samples)
    sampling = random.sample(sampling_list, len(sampling_list))
    genotype = genotype[sampling, :]
    new_samples = list()
    for i in range(len(sampling)):
        new_samples.append(samples[sampling[i]])
    samples = new_samples
    phenotype = list()
    for sample in samples:
        phenotype.append(phenotype_lib[sample])
    population = list()
    for sample in samples:
        population.append(population_lib[sample])
    strip_index_list = load_strip_list(num_samples)
    g_trainings, g_tests = strip_samples(genotype, strip_index_list)
    p_trainings, p_tests = strip_samples(phenotype, strip_index_list)
    pop_trainings, pop_tests = strip_samples(population, strip_index_list)
    for CV in range(CV_fold):
        g_training = g_trainings[CV]; g_test = g_tests[CV]
        p_training = p_trainings[CV]; p_test = p_tests[CV]
        pop_test = pop_tests[CV]; 
        LinearRegression.linear_regression(g_training, g_test, p_training, 
                                           p_test, pop_test, CV)


def get_genotype():
    lines = open(genotype_file, "r").readlines()
    num_samples = len(lines[0][:-2].split("\t"))-5
    num_markers = len(lines)-1
    genotype = np.zeros([num_samples, num_markers])
    samples = list(); markers = list()
    for i,line in enumerate(lines):
        temp_list = line[:-2].split("\t")
        if(i == 0):
            for temp in temp_list[5:]:
                samples.append(str(temp))
        else:
            markers.append(temp_list[:5])
            for j,temp in enumerate(temp_list[5:]):
                genotype[j,i-1] = float(temp)
    return markers, samples, genotype

def get_phenotype(samples):
    lines = open(phenotype_file, "r").readlines()
    phenotype_lib = {}
    for line in lines[1:]:
        temp_list = line[:-1].split(" ")
        temp_name = str(temp_list[0])
        temp_val = float(temp_list[2])
        for sample in samples:
            if(temp_name+"@" in sample):
                phenotype_lib.update({sample:temp_val})
    return phenotype_lib
        
def get_population(samples):
    lines = open("data/population.txt", "r").readlines()
    IGRC_list = list(); NSFTV_list = list(); population_list = list()
    for line in lines:
        temp_list = line.split("\t")
        IGRC_list.append(temp_list[4])
        NSFTV_list.append(temp_list[5])
        population_list.append(temp_list[8])
    population_lib = {}
    for sample in samples:
        for IGRC, NSFTV, population in zip(IGRC_list, NSFTV_list,
                population_list):
            if(IGRC + "@" in sample or NSFTV + "@" in sample):
                population_lib.update({sample:population})
    return population_lib


def load_strip_list(num_samples):
    split_size = int(num_samples/CV_fold)
    strip_index_list = list()
    for i in range(1, CV_fold):
        strip_index_list.append(split_size * i)
    strip_index_list.append(num_samples)
    return strip_index_list


def strip_samples(genotype, strip_index_list):
    trainings = list(); tests = list()
    for CV in range(CV_fold):
        if(CV == 0):
            start = 0
        else:
            start = strip_index_list[CV-1]
        end = strip_index_list[CV]
        test = genotype[start:end]
        tests.append(test)
        training = np.r_[genotype[:start],genotype[end:]]
        trainings.append(training)
    return trainings, tests

if __name__ == "__main__":
    genotype_file = "data/genotype_ppca.tsv"
    phenotype_file = "data/phenotype.txt"
    CV_fold = 3
    main()
