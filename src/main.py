#!/usr/bin/env python3
import sys
import os

# Add the current directory to the Python path to find the module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.data_pipeline.main import main_data_pipeline
from src.knowledge_graph_builder.main import main_kg_builder


def main(query, output_folder, config_path="./config.toml"):
    """
    Main function to run the data pipeline and knowledge graph builder.

    Args:
        query (str): The query string to be processed.
        output_folder (str): The folder where the output files will be saved.
    """
    experiment_name = "pubmed_" + query.replace(" ", "_").lower()

    # Run the data pipeline
    results = main_data_pipeline(
        user_query=query, experiment_name=experiment_name, output_folder=output_folder
    )

    pdf_download_folder = results.get("pdf_download_folder", None)

    # Build the knowledge graph
    main_kg_builder(
        input_path=pdf_download_folder,
        output_dir=output_folder,
        config_path=config_path,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the data pipeline and knowledge graph builder."
    )
    parser.add_argument("--query", type=str, help="The query string to be processed.")
    parser.add_argument(
        "--output_folder",
        type=str,
        default="/home/exouser/masters-thesis/ai-knowledge-graph-main/data2",
        help="The folder where the output files will be saved.",
    )
    parser.add_argument(
        "--config_path",
        type=str,
        default="./config.toml",
        help="Path to the configuration file.",
    )

    args = parser.parse_args()

    main(args.query, args.output_folder, args.config_path)
