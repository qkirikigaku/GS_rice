import numpy as np
import math

def do_ppca(genotype, n_components):
    missed_flag = find_missed(genotype)
    num_sample = len(genotype)
    num_marker = len(genotype[0])
    Z, W, Mu, sigma = initialize(num_sample, num_marker, n_components)
    LLs = list()
    old_LL = -1e10
    for iter in range(100):
        print("iter:",iter)
        genotype = fill_missed(genotype, missed_flag, Z, W, Mu)
        Z, ZZ = E_step(genotype, Z, W, Mu, sigma)
        Mu, W, sigma = M_step(genotype, Z, ZZ, W, Mu, sigma)
        temp_LL = calc_likelihood(genotype, Z, ZZ, W, Mu, sigma)
        LLs.append(temp_LL)
        improved = temp_LL - old_LL
        print("  ", improved, "improved")
        if(improved < 1e-4):
            break
        else:
            old_LL = temp_LL
    genotype = fill_missed(genotype, missed_flag, Z, W, Mu)
    write_data(genotype, Z, W, Mu, sigma, LLs)

def find_missed(genotype):
    num_sample = len(genotype); num_marker = len(genotype[0])
    missed_flag = np.zeros([num_sample, num_marker])
    for i in range(num_sample):
        for j in range(num_marker):
            if(np.isnan(genotype[i,j])):
                missed_flag[i,j] = 1
    return missed_flag

def initialize(num_sample, num_marker, n_components):
    Z = np.random.normal(0, 1, (num_sample, n_components))
    W = np.random.normal(0, 0.1, (num_marker, n_components))
    Mu = np.random.normal(0, 1, num_marker)
    sigma = np.random.uniform(0, 1)
    return Z, W, Mu, sigma

def fill_missed(genotype, missed_flag, Z, W, Mu):
    num_sample = len(genotype); num_marker = len(genotype[0])
    for i in range(num_sample):
        Z_i = Z[i]
        for j in range(num_marker):
            W_j = W[j]
            Mu_j = Mu[j]
            if(missed_flag[i,j] == 1):
                genotype[i,j] = Z_i @ W_j + Mu_j
    return genotype

def E_step(genotype, Z, W, Mu, sigma):
    num_sample = len(genotype); num_marker = len(genotype[0])
    n_components = len(W[0])
    ZZ = np.zeros([num_sample, n_components, n_components])
    M = W.T @ W + np.diag([sigma]*n_components)
    for i in range(num_sample):
        Z[i] = np.linalg.inv(M) @ W.T @ (genotype[i] - Mu)
        ZZ[i] = sigma * np.linalg.inv(M) + np.outer(Z[i], Z[i])
    return Z, ZZ

def M_step(genotype, Z, ZZ, W, Mu, sigma):
    num_sample = len(genotype); num_marker = len(genotype[0])
    n_components = len(W[0])
    Mu = np.zeros([num_marker])
    for i in range(num_sample):
        Mu += genotype[i]/num_sample
    first_component = np.zeros([num_marker, n_components])
    second_component = np.zeros([n_components, n_components])
    for i in range(num_sample):
        first_component += np.outer((genotype[i] - Mu), Z[i])
        second_component += ZZ[i]
    W = first_component @ np.linalg.inv(second_component)
    sigma = 0
    for i in range(num_sample):
        sigma += np.dot((genotype[i] - Mu), (genotype[i] - Mu)) -\
                 2 * Z[i].T @ W.T @ (genotype[i] - Mu) +\
                 np.trace(ZZ[i] @ W.T @ W)
    sigma /= num_sample * num_marker
    return Mu, W, sigma

def calc_likelihood(genotype, Z, ZZ, W, Mu, sigma):
    temp_LL = 0
    num_sample = len(genotype); num_marker = len(genotype[0])
    n_components = len(W[0])
    for i in range(num_sample):
        temp_LL -= (num_marker/2) * np.log(2 * math.pi * sigma) +\
                   (1/2) * np.trace(ZZ[i]) +\
                   (1/(2*sigma)) * np.dot((genotype[i]-Mu), (genotype[i]-Mu)) -\
                   (1/sigma) * Z[i].T @ W.T @ (genotype[i] - Mu) +\
                   (1/(2*sigma)) * np.trace(ZZ[i] @ W.T @ W) +\
                   (n_components/2) * np.log(2 * math.pi)
    return temp_LL

def write_data(genotype, Z, W, Mu, sigma, LLs):
    num_sample = len(genotype); num_marker = len(genotype[0])
    n_components = len(W[0])
    output = open("result/ppca/genotype.tsv", "w")
    for i in range(num_sample):
        for j in range(num_marker):
            output.write(str(genotype[i,j]) + "\t")
        output.write("\n")
    output = open("result/ppca/Z.tsv", "w")
    for i in range(num_sample):
        for j in range(n_components):
            output.write(str(Z[i,j]) + "\t")
        output.write("\n")
    output = open("result/ppca/W.tsv", "w")
    for i in range(num_marker):
        for j in range(n_components):
            output.write(str(W[i,j]) + "\t")
        output.write("\n")
    output = open("result/ppca/MuSigma.txt", "w")
    for i in range(num_marker):
        output.write(str(Mu[i]) + "\t")
    output.write("\n")
    output.write(str(sigma) + "\n")
    output = open("result/ppca/LLs.txt", "w")
    for i in range(len(LLs)):
        output.write(str(LLs[i]) + "\n")


