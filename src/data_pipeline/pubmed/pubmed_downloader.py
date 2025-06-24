#!/usr/bin/env python
"""
PubMed Central PDF Fetcher

This tool downloads PDF articles from various publishers using PubMed IDs.
It supports multiple publishers including NEJM, Science Direct, and PubMed Central.
"""

import os
import sys
import argparse
import requests
import pandas as pd
from typing import Optional, Tuple
from Bio import Entrez
import xml.etree.ElementTree as ET
import time
from dotenv import load_dotenv
from indra.literature.pmc_client import extract_text

# Load environment variables from .env file
load_dotenv()

# Set your email for Entrez from the environment variable
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL")
if not ENTREZ_EMAIL:
    raise ValueError("Please set ENTREZ_EMAIL in your .env file.")
Entrez.email = ENTREZ_EMAIL

# Optionally set the PubMed API key if provided
PUBMED_API_KEY = os.getenv("PUBMED_API_KEY")
if PUBMED_API_KEY:
    print("Using PubMed API key for enhanced rate limits.")
    Entrez.api_key = PUBMED_API_KEY


def check_full_text_availability(pmid: str) -> Tuple[bool, Optional[str]]:
    """Check if full text is available in PMC and get PMC ID if it exists.

    Args:
        pmid: PubMed ID of the article

    Returns:
        Tuple of (availability boolean, PMC ID if available)
    """
    try:
        print(f"Checking PMC availability for PMID {pmid}")
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)

        if not handle:
            print(f"No PMC link found for PMID {pmid}")
            return False, None

        xml_content = handle.read()
        handle.close()

        # Parse XML to get PMC ID
        root = ET.fromstring(xml_content)
        linksetdb = root.find(".//LinkSetDb")
        if linksetdb is None:
            print(f"No PMC ID found for PMID {pmid}")
            return False, None

        id_elem = linksetdb.find(".//Id")
        if id_elem is None:
            print(f"No PMC ID element found for PMID {pmid}")
            return False, None

        pmc_id = id_elem.text
        print(f"Found PMC ID {pmc_id} for PMID {pmid}")
        return True, pmc_id

    except Exception as e:
        print(f"Error checking PMC availability for PMID {pmid}: {str(e)}")
        return False, None


def get_full_text(pmid: str) -> Optional[str]:
    """Get full text of the article if available through PMC.

    Handles truncated responses by making additional requests.

    Args:
        pmid: PubMed ID of the article

    Returns:
        Full text content if available, None otherwise
    """
    try:
        # First check availability and get PMC ID
        available, pmc_id = check_full_text_availability(pmid)
        if not available or pmc_id is None:
            print(f"Full text not available in PMC for PMID {pmid}")
            return None

        print(f"Fetching full text for PMC ID {pmc_id}")
        content = ""
        retstart = 0

        while True:
            full_text_handle = Entrez.efetch(
                db="pmc", id=pmc_id, rettype="xml", retstart=retstart
            )

            if not full_text_handle:
                break

            chunk = full_text_handle.read()
            full_text_handle.close()

            if isinstance(chunk, bytes):
                chunk = chunk.decode("utf-8")

            content += chunk

            # Check if there might be more content
            if "[truncated]" not in chunk and "Result too long" not in chunk:
                break

            # Increment retstart for next chunk
            retstart += len(chunk)

            # Add small delay to respect API rate limits
            time.sleep(0.5)

        result = extract_text(content)
        return result

    except Exception as e:
        print(f"Error getting full text for PMID {pmid}: {str(e)}")
        return None


def parse_arguments():
    """Parse and validate command line arguments."""
    parser = argparse.ArgumentParser()
    parser._optionals.title = "Flag Arguments"
    parser.add_argument(
        "-pmids",
        help="Comma separated list of pmids to fetch. Must include -pmids or -pmf.",
        default="%#$",
    )
    parser.add_argument(
        "-csv",
        help="File with pmids to fetch inside, one pmid per line. Optionally, the file can be a tsv with a second column of names to save each pmid's article with (without '.pdf' at the end). Must include -pmids or -pmf",
        default="%#$",
    )
    parser.add_argument(
        "-out",
        help="Output directory for fetched articles. Default: fetched_pdfs",
        default="fetched_pdfs",
    )
    parser.add_argument(
        "-errors",
        help="Output file path for pmids which failed to fetch. Default: unfetched_pmids.tsv",
        default="unfetched_pmids.tsv",
    )

    args = vars(parser.parse_args())

    # Check for required arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args["pmids"] == "%#$" and args["csv"] == "%#$":
        print("Error: Either -pmids or -pmf must be used. Exiting.")
        sys.exit(1)

    if args["pmids"] != "%#$" and args["csv"] != "%#$":
        print("Error: -pmids and -pmf cannot be used together. Ignoring -pmf argument")
        args["pmf"] = "%#$"

    return args


def run_pubmed_downloader(output_dir, input_csv_file, errors_csv_file):
    """Main execution function."""
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        print(f"Output directory of {output_dir} did not exist. Created the directory.")
        os.mkdir(output_dir)

    df = pd.read_csv(input_csv_file)
    pmids = df["PMID"].tolist()
    names = df["Title"].tolist()
    names = [str(name) for name in names]  # Ensure all names are strings

    names = [
        name.replace(" ", "_").replace(" ", "_").replace(".", "_").replace("/", "_")
        for name in names
    ]

    # Initialize errors CSV file with header if it doesn't exist
    if not os.path.exists(errors_csv_file):
        with open(errors_csv_file, 'w') as f:
            f.write("pmid,name\n")

    # Process each PMID
    for pmid, name in zip(pmids, names):
        print(f"Trying to fetch pmid {pmid}")

        try:
            output = get_full_text(str(pmid))
            if output:
                # Save the full text to a file
                filepath = os.path.join(output_dir, f"{name}.txt")
                with open(filepath, "w", encoding="utf-8") as f:
                    f.write(output)
                print(f"** fetching of reprint {pmid} succeeded")
            else:
                print(f"** fetching of reprint {pmid} failed, no full text available")
                # Write error immediately to CSV
                with open(errors_csv_file, 'a') as f:
                    f.write(f"{pmid},{name}\n")
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                print(f"** fetching of reprint {pmid} failed, article not found")
                # Write error immediately to CSV
                with open(errors_csv_file, 'a') as f:
                    f.write(f"{pmid},{name}\n")
        except Exception as e:
            print(f"** fetching of reprint {pmid} failed with error: {str(e)}")
            # Write error immediately to CSV
            with open(errors_csv_file, 'a') as f:
                f.write(f"{pmid},{name}\n")


if __name__ == "__main__":
    args = parse_arguments()
    run_pubmed_downloader(args["out"], args["csv"], args["errors"])
