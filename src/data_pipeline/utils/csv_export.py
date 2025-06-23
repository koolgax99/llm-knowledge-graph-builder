"""
xlsx_export.py: Contains functions to export the search results and metadata to an Excel file with two tabs.
"""

import os
import pandas as pd
from src.data_pipeline.utils.database import SearchResult


def check_filename(filename):
    """
    Check if the filename is valid.
    - Filename should not be empty.
    - Filename should not contain slashes.
    - Filename should end with .xlsx
    - If no extension is provided, append .xlsx.
    - If no folder path is provided, use "data" as the default folder.
    - Ensure the output folder exists.
    - If filename does not contain folder path, prepend it.
    - Return the full path of the filename.
    """
    if not filename:
        raise ValueError("Filename cannot be empty.")

    if not filename.endswith(".csv"):
        filename += ".csv"

    if not os.path.splitext(filename)[1]:
        filename += ".csv"

    # check if folder exists for filename
    folder = os.path.dirname(filename)
    if not os.path.exists(folder):
        os.makedirs(folder)

    return filename


def export_to_csv(session, metadata_id, search_filename="output.csv"):
    """
    Export search results and metadata to an Excel file with two sheets.
    Sheet1: Search Results
    Sheet2: Metadata
    """
    # Retrieve data from database
    results = (
        session.query(SearchResult)
        .filter(SearchResult.metadata_id == metadata_id)
        .all()
    )

    # Create DataFrame for search results
    results_data = []
    for result in results:
        results_data.append(
            {
                "PMID": result.pmid,
                "Title": result.title,
                "Authors": result.authors,
                "Abstract": result.abstract,
                "DOI": result.doi,
                "Link": result.link,
                "Year": result.year,
            }
        )
    df_results = pd.DataFrame(results_data)

    # Check the filename
    search_filename = check_filename(search_filename)

    # Save each dataframe as a separate CSV file
    df_results.to_csv(search_filename, index=False)
