
libary_file = open("libray.tsv", "r")


def prepare_libary():
    LIBARY = {}
    for line in libary_file.readlines():
        line = line.split(",")
        if line[1] != "cell type":
            LIBARY[line[0]] = line[1].strip("\n")
    libary_file.close()
    return LIBARY

def check_gene(gene: str, LIBARY):
    if gene in LIBARY:
        return LIBARY[gene]
    else:
        return "Gene not found in libary"

def validate_cluster_genes(cluster_genes: list, LIBARY):   
    cell_types = []
    for item in cluster_genes:
        cell_types.append(check_gene(item, LIBARY))

    count_dict = {}
    for item in cell_types:
        if item in count_dict:
            count_dict[item] += 1
        else:
            count_dict[item] = 1
    return count_dict

def read_cluster_genes(cluster_file, cluster_def_name_col, gene_def_name_col):
    cluster_file = open(cluster_file, "r")
    parameters  = cluster_file.readline().split("\t")
    index_cluster = parameters.index(cluster_def_name_col)
    index_gene = parameters.index(gene_def_name_col)
    cluster_gene = {}
    for line in cluster_file.readlines():
        line = line.split("\t")
        cluster = line[index_cluster]
        gene = line[index_gene]
        if cluster in cluster_gene:
            cluster_gene[cluster].append(gene)
        else:
            cluster_gene[cluster] = [gene]
    cluster_file.close()
    return cluster_gene # change to list of tuple instad of dict 

def annotate_clusters(cluster_gene_dict, LIBARY):
    annotated_clusters = {}
    for cluster in cluster_gene_dict:
        annotated_clusters[cluster] = validate_cluster_genes(cluster_gene_dict[cluster], LIBARY)
    return annotated_clusters

LIBARY = prepare_libary()
gene_def_name_col = "seurat_clusters"
cluster_def_name_col = "markers"

cluster_gene_dict = read_cluster_genes("cluster_genes.txt", cluster_def_name_col, gene_def_name_col)
annotate_clusters(cluster_gene_dict, LIBARY)
