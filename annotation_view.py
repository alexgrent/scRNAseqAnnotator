import annotation_libary as annoation

LIBARY = annoation.prepare_libary()


gene_def_name_col = input("Column of gene: ")
cluster_def_name_col = input("Column of cluster: ")
file_dir = input("File directory (top ten marker file): ")

cluster_gene_dict = annoation.read_cluster_genes(file_dir, cluster_def_name_col, gene_def_name_col)
annoation = annoation.annotate_clusters(cluster_gene_dict, LIBARY)
print(annoation)