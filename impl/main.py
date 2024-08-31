import json
from flask import Flask, request, jsonify
from services import SCRNAProcessService  # Assuming this is your service module

app = Flask(__name__)


@app.route("/", methods=['GET'])
def init_application():
    return "Application starting"


@app.route("/cluster", methods=['POST'])
def get_clustering_per_resolution():
    data = request.get_json()

    if 'resolution' not in data:
        return jsonify({"error": "Resolution parameter missing"}), 400

    # Assuming SCRNAProcessService.load_clustering_per_resolution handles clustering
    SCRNAProcessService.load_clustering_per_resolution(data['resolution'])

    return jsonify({"message": "Clustering initiated for resolution: {}".format(data['resolution'])}), 200


@app.route("/process", methods=['POST'])
def process_scDATA():
    data = request.get_json()

    required_keys = ["scObjectFilePath", "min_genes", "min_cells", "top_genes", "batch_key", "pc", "resolution"]
    if not all(key in data for key in required_keys):
        return jsonify({"error": "Missing one or more required parameters"}), 400

    # Assuming SCRNAProcessService.process_scRNAdata handles the scRNA processing
    SCRNAProcessService.process_scRNAdata(
        data["scObjectFilePath"],
        data["min_genes"],
        data["min_cells"],
        data["top_genes"],
        data["batch_key"],
        data["pc"],
        data["resolution"]
    )

    return jsonify({"status": "processing"}), 200



if __name__ == "__main__":
    app.run(debug=True)
