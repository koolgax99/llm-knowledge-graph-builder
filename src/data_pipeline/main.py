"""
main.py: Core processing module for the AI agent that retrieves PubMed articles based on a user query.
"""

import logging
import sys
from src.data_pipeline.utils.extract_values import extract_years_from_query
from src.data_pipeline.pubmed.pubmed_search import run_pubmed_search
from src.data_pipeline.pubmed.pubmed_downloader import run_pubmed_downloader
from src.data_pipeline.utils.database import (
    init_db,
    store_metadata,
    store_search_results,
    get_engine_session,
)
from src.data_pipeline.utils.csv_export import export_to_csv
from src.data_pipeline.utils.constants import DEFAULT_DATE_RANGE
from dotenv import load_dotenv

# from utils.pdf_to_text import main_pdf_converter
import os

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)


def main_data_pipeline(
    user_query=None,
    min_year=None,
    max_year=None,
    experiment_name=None,
    output_folder=None,
    max_results=None
):
    """
    Process a PubMed search request with the given parameters.

    Args:
        user_query (str, optional): The search query text. If None, will use a default prompt.
        min_year (int, optional): Minimum publication year. If None, will try to extract from query or use default.
        max_year (int, optional): Maximum publication year. If None, will try to extract from query or use default.
        export_results (bool, optional): Whether to export results to CSV. Defaults to False.

    Returns:
        dict: A dictionary containing the search results and metadata
    """
    # Get user query if not provided
    if user_query is None:
        user_query = input("Enter your search query: ")

    user_query = user_query.strip()

    # Extract date range from query if not provided
    if min_year is None or max_year is None:
        extracted_years = extract_years_from_query(user_query)
        if extracted_years:
            min_year, max_year, user_query = extracted_years
            logging.info(
                "Extracted date range from query: %d to %d", min_year, max_year
            )
        else:
            if max_year is None:
                max_year = DEFAULT_DATE_RANGE["MAX_YEAR"]
            if min_year is None:
                min_year = DEFAULT_DATE_RANGE["MIN_YEAR"]
            logging.info("Using date range: %d to %d", min_year, max_year)

    # Execute PubMed search
    logging.info("Executing PubMed search...")
    search_results, final_query = run_pubmed_search(
        user_query, min_year, max_year, ret_max=max_results
    )
    logging.info("Retrieved %d search results", len(search_results))

    # Initialize database and store results
    engine = init_db()
    session = get_engine_session(engine)
    metadata_id = None

    try:
        metadata_id = store_metadata(session, min_year, max_year, final_query)
        store_search_results(session, search_results, metadata_id)
        session.commit()
        logging.info("Data stored successfully in the database.")
    except Exception as e:
        session.rollback()
        logging.error("Error storing data: %s", str(e))
    finally:
        session.close()

    # Export to CSV if requested
    if metadata_id:
        try:
            logging.info("Exporting results to CSV...")
            pmid_csv_filename = os.path.join(
                output_folder or "data", f"{experiment_name}_pubmed_results.csv"
            )
            logging.info("Exporting results to %s", pmid_csv_filename)
            export_to_csv(session, metadata_id, pmid_csv_filename)
            logging.info("Data exported successfully to %s", pmid_csv_filename)
        except Exception as e:
            logging.error("Error exporting data: %s", str(e))

    # Call the Pubmed Downloader to download full articles
    try:
        logging.info("Running PubMed downloader...")
        pdf_download_folder = os.path.join(
            output_folder or "data", f"{experiment_name}_pubmed_downloads"
        )
        os.makedirs(pdf_download_folder, exist_ok=True)
        errors_csv_filename = os.path.join(
            output_folder or "data", f"{experiment_name}_pubmed_download_errors.csv"
        )
        run_pubmed_downloader(
            pdf_download_folder, pmid_csv_filename, errors_csv_filename
        )
    except Exception as e:
        logging.error("Error during PubMed downloader: %s", str(e))

    # Convert the pdf files to text
    # try:
    #     logging.info("Converting PDF files to text...")
    #     txt_converted_folder = os.path.join(output_folder or "data", f"{experiment_name}_txt_converted")
    #     main_pdf_converter(pdf_download_folder, txt_converted_folder or "data")
    #     logging.info("PDF files converted to text and saved in %s", txt_converted_folder)
    # except Exception as e:
    #     logging.error("Error converting PDF files to text: %s", str(e))

    # Return results for potential further processing
    return {
        "search_query": user_query,
        "date_range": {"min_year": min_year, "max_year": max_year},
        "mesh_strategy": final_query,
        "results_count": len(search_results),
        "metadata_id": metadata_id,
        "output_filename": pmid_csv_filename,
        "errors_filename": errors_csv_filename,
        "pdf_download_folder": pdf_download_folder,
        # "txt_converted_folder": txt_converted_folder
    }


# This is maintained for backward compatibility when running as a script
if __name__ == "__main__":
    import argparse

    # Command line interface for direct invocation
    parser = argparse.ArgumentParser(description="AI agent for PubMed searches")
    parser.add_argument("--query", type=str, help="User query for the search")
    parser.add_argument("--min_year", type=int, help="Minimum publication year")
    parser.add_argument("--max_year", type=int, help="Maximum publication year")
    parser.add_argument(
        "--experiment_name",
        type=str,
        default="pubmed_search",
        help="Name of the experiment for output files",
    )
    parser.add_argument(
        "--output_folder",
        type=str,
        default="/home/exouser/masters-thesis/ai-knowledge-graph-main/test_data",
        help="Folder to save output files",
    )
    args = parser.parse_args()

    try:
        main_data_pipeline(
            user_query=args.query,
            min_year=args.min_year,
            max_year=args.max_year,
            experiment_name=args.experiment_name,
            output_folder=args.output_folder,
        )
    except KeyboardInterrupt:
        sys.exit("Process interrupted by user.")
