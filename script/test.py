def main():
    count = 0
    lines = open(genotype_file, "r").readlines()
    for line in lines[1:]:
        temp_list = line.split("\t")
        for temp in temp_list[5:]:
            if(temp == "N"):
                count += 1
    print("count: ",count)

if __name__ == "__main__":
    genotype_file = "data/genotype_reshaped.tsv"
    phenotype_file = "data/phenotype.txt"
    main()
