def main():
    columns,markers = get_columns()
    genotype = get_genotype()
    output_file = "data/genotype_ppca.tsv"
    output = open(output_file, "w")
    for i in range(len(markers)+1):
        if(i == 0):
            for column in columns:
                output.write(str(column) + "\t")
        else:
            for marker in markers[i-1]:
                output.write(str(marker) + "\t")
            for gene in genotype:
                output.write(str(gene[i-1]) + "\t")
        output.write("\n")


def get_columns():
    columns = open(genotype_columns_file, "r").readline().split("\t")
    columns[len(columns)-1] = columns[len(columns)-1][:-1]
    lines = open(genotype_columns_file, "r").readlines()
    markers = list()
    for line in lines[1:]:
        temp = line.split("\t")[0:5]
        markers.append(temp)
    return columns, markers

def get_genotype():
    genotype = list()
    lines = open(genotype_file, "r").readlines()
    for i,line in enumerate(lines):
        genotype.append([])
        temp_list = line.split("\t")
        for temp in temp_list[:-1]:
            genotype[i].append(float(temp))
    return genotype

if __name__ == "__main__":
    genotype_file = "result/ppca/genotype.tsv"
    genotype_columns_file = "data/genotype_reshaped.tsv"
    main()
