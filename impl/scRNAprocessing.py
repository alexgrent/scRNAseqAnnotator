import scanpy as sc
import anndata as ad

import pooch

class scRNAAnalysis:

    

    def __init__(self, scObject, min_genes = 100, min_cells = 3, top_genes = 2000, batch_key = "sample", pc = 50, resolution = 0.5):
        self.leiden_clustering = None
        self.umpa_plot = None
        self.pc_plot = None
        self.qc_scatter = None
        self.qc_violin = None
        self.scObject = scObject
        self.min_genes = min_genes
        self.min_cells = min_cells
        self.top_genes = top_genes
        self.batch_key = batch_key
        self.pc = pc
        self.res = resolution

    def run_analysis(self):
        # mitochondrial genes, "MT-" for human, "Mt-" for mouse
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        # ribosomal genes
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        # hemoglobin genes
        adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
        )

        self.qc_violin = sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
        )

        self.qc_scatter = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

        # QC with cell filtering, gene filtering and doublet detection
        sc.pp.filter_cells(adata, min_genes=self.min_genes)
        sc.pp.filter_genes(adata, min_cells=self.min_cells)
        #sc.pp.scrublet(adata, batch_key="sample")

        print("Normalization")
        # Normalization
        adata.layers["counts"] = adata.X.copy()
        # Normalizing to median total counts
        sc.pp.normalize_total(adata)
        # Logarithmize the data
        sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(adata, n_top_genes=self.top_genes, batch_key=self.batch_key)
        sc.tl.pca(adata)
        sc.pl.pca_variance_ratio(adata, n_pcs=self.pc, log=True)

        print("PCA")
        # show PCA plots
        self.pc_plot = sc.pl.pca(
            adata,
            color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
            dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
            ncols=2,
            size=2,
        )

        # clustering and UMAP
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

        self.umpa_plot = sc.pl.umap(
            adata,
            color="sample",
            # Setting a smaller point size to get prevent overlap
            size=2,
        )

        print("clustering")
        # Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
        sc.tl.leiden(adata, key_added=f"leiden_res_{self.res:4.2f}", resolution= self.res, flavor="igraph")
        self.leiden_clustering = sc.pl.umap(adata, color=[f"leiden_res_{self.res:4.2f}"])





EXAMPLE_DATA = pooch.create(
    path=pooch.os_cache("scverse_tutorials"),
    base_url="doi:10.6084/m9.figshare.22716739.v1/",
)
EXAMPLE_DATA.load_registry_from_doi()

samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
}
adatas = {}

for sample_id, filename in samples.items():
    path = EXAMPLE_DATA.fetch(filename)
    sample_adata = sc.read_10x_h5(path)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
print(adata.obs["sample"].value_counts())

run = scRNAAnalysis(adata, 100,3, 2000, "sample",  50, 0.5)
run.run_analysis()
print(run)
