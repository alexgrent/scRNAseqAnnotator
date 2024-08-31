from scRNAprocessing import scRNAAnalysis

class SCRNAProcessService:


    @staticmethod
    def process_scRNAdata(scFilePath, min_genes, min_cells, top_genes, batch_key, pc, resolution):
        analyse_process = scRNAAnalysis(scFilePath, min_genes, min_cells, top_genes, batch_key, pc, resolution)


    @staticmethod
    def load_clustering_per_resolution(resolution: float):
         return 4