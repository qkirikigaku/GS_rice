import pandas as pd
import numpy as np
import time

non_defect_rate = 0.99
MAF_threshold = 0.05

valid_base_list = ["R","M","W","S","Y","K","A","C","G","T"]
hetero_list = ["R", "M", "W", "S", "Y", "K"]

def main():
    genotype = pd.read_csv(genotype_file, sep="\t")
    drop_list = list()
    genotype.reset_index(drop=True, inplace=True)
    columns_list = list(genotype.columns)
    genotype = genotype.values
    genotype = filter_minor_allele_exist(genotype)
    genotype = pd.DataFrame(data=genotype, columns=columns_list)
    genotype.reset_index(drop=True, inplace=True)
    genotype.to_csv("data/genotype_shaped.tsv", sep="\t", index=False)


def filter_minor_allele_exist(genotype):
    drop_list = list()
    for index in range(len(genotype)):
        temp_allele = list(genotype[index, 1])
        if(len(temp_allele) != 3):
            drop_list.append(index)
        else:
            base_candidate = list()
            base_candidate.append(temp_allele[0])
            base_candidate.append(temp_allele[2])
            bases = list(genotype[index, 11:])
            valid_count = 0
            for base in bases:
                if(base in valid_base_list):
                    valid_count += 1
            ratio = float(valid_count/len(bases))
            base_dic = {}
            for base in bases:
                if(base not in base_dic and base in valid_base_list):
                    base_dic.update({base:1})
                elif(base in valid_base_list):
                    val = base_dic[base]
                    val += 1
                    base_dic.update({base:val})
            base_count = [0,0]
            for key in base_dic.keys():
                if(key in hetero_list):
                    base_count[0] += base_dic[key]
                    base_count[1] += base_dic[key]
                elif(key == base_candidate[0]):
                    base_count[0] += base_dic[key]
                elif(key == base_candidate[1]):
                    base_count[1] += base_dic[key]
            MAF = min(base_count[0], base_count[1])
            minor_allele = base_candidate[np.argmin(base_count)]
            MAF /= (2*valid_count)
            if(ratio <= non_defect_rate or MAF <= MAF_threshold):
                drop_list.append(index)
            else:
                temp_bases = list()
                for temp_base in bases:
                    new_base = transform_base_to_num(minor_allele, temp_base)
                    temp_bases.append(new_base)
                genotype[index, 11:] = temp_bases
        if(index % 700 == 0): print(str(float(index/len(genotype))*100) + "%")
    genotype = np.delete(genotype, drop_list, 0)
    return genotype


def transform_base_to_num(minor_allele, base):
    if(base == minor_allele): return -1.0
    elif(base in hetero_list): return 0.0
    elif(base != "N"): return 1.0
    else: return "N"

if __name__ == "__main__":
    genotype_file = "data/genotype.txt"
    phenotype_file = "data/phenotype.txt"
    main()
