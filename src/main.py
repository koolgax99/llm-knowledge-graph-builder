#!/usr/bin/env python3
import sys
import os

# Add the current directory to the Python path to find the module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.data_pipeline.main import main_data_pipeline
from src.knowledge_graph_builder.main import main_kg_builder


def main(query, max_results, output_folder=None, config_path="./config.toml"):
    """
    Main function to run the data pipeline and knowledge graph builder.

    Args:
        query (str): The query string to be processed.
        output_folder (str): The folder where the output files will be saved.
    """
    if output_folder is None:
        output_folder = os.path.abspath(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output")
        )

    experiment_name = "pubmed_" + query.replace(" ", "_").lower()
    experiment_name = experiment_name.replace("(", "").replace(")", "").replace(",", "_")

    # Run the data pipeline
    results = main_data_pipeline(
        user_query=query,
        experiment_name=experiment_name,
        output_folder=output_folder,
        max_results=max_results,
    )

    #pdf_download_folder = results.get("pdf_download_folder", None)

    pdf_download_folder = os.path.join(output_folder, f"{experiment_name}_pubmed_downloads")

    # Build the knowledge graph
    main_kg_builder(
        input_path=pdf_download_folder,
        output_dir=output_folder,
        config_path=config_path,
        experiment_name=experiment_name
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the data pipeline and knowledge graph builder."
    )
    parser.add_argument("--query", type=str, help="The query string to be processed.")
    parser.add_argument(
        "--max_results",
        type=int,
        default=20,
        help="The maximum number of results to return.",
    )
    parser.add_argument(
        "--output_folder",
        type=str,
        default=None,
        help="The folder where the output files will be saved.",
    )
    parser.add_argument(
        "--config_path",
        type=str,
        default="./config.toml",
        help="Path to the configuration file.",
    )

    args = parser.parse_args()

    main(args.query, args.max_results, args.output_folder, args.config_path)
