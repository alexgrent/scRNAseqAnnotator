


file_pangloa = open("PanglaoDB_markers_27_Mar_2020.tsv", "r")
libray = open("libray.tsv", "w")

for line in file_pangloa.readlines():
    line = line.split("\t")
    print(line)
    gene_symbol = line[1]
    cell_type = line[2]
    print(gene_symbol + "," + cell_type + "\n")
    libray.write(gene_symbol + "," + cell_type + "\n")

file_pangloa.close()
libray.close()
