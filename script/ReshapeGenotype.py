import pandas as pd


def main():
    genotype = pd.read_csv(genotype_file, sep="\t")
    drop_list = ["assembly#", "center", "protLSID", "assayLSID",\
                 "panelLSID", "QCcode"]
    genotype.drop(drop_list, axis=1, inplace=True)
    sample_candidates = list(genotype.columns[5:])
    new_samples = list()
    lines = open(phenotype_file, "r").readlines()
    for line in lines[1:]:
        temp_list = line.split(" ")
        temp_sample = str(temp_list[0])
        for sample_candidate in sample_candidates:
            if(temp_sample+"@" in str(sample_candidate)):
                new_samples.append(sample_candidate)
                break
    print("sample_candidate: " + str(len(sample_candidates)) + " samples")
    print("new_sample: " + str(len(new_samples)) + " samples")
    drop_samples = list(set(sample_candidates).difference(set(new_samples)))
    print("drop_sample: " + str(len(drop_samples)) + " samples")
    genotype.drop(drop_samples, axis=1, inplace=True)
    genotype.to_csv("data/genotype_reshaped.tsv", sep="\t", index=False)

if __name__ == "__main__":
    genotype_file = "data/genotype_shaped.tsv"
    phenotype_file = "data/phenotype.txt"
    main()
