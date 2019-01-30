import numpy as np
from sklearn.decomposition import PCA
import PPCA

def main():
    lines = open(genotype_file, "r").readlines()
    columns = lines[0].split("\t")
    sample_names = columns[5:]
    num_sample = len(sample_names)
    num_feature = len(lines) - 1
    genotype = np.zeros([num_sample, num_feature])
    for i,line in enumerate(lines[1:]):
        temp_list = line.split("\t")[5:]
        for j,temp in enumerate(temp_list):
            if(temp in ["N", "N\n"]):
                genotype[j,i] = "NaN"
            else:
                genotype[j,i] = float(temp)
    PPCA.do_ppca(genotype, 8)

if __name__ == "__main__":
    genotype_file = "data/genotype_reshaped.tsv"
    main()
