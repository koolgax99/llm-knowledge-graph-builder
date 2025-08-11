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
import re
import concurrent.futures
import threading

# Load environment variables from .env file
load_dotenv()

lock = threading.Lock()

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
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", link_name="pubmed_pmc", id=pmid, retmode='text')

        handle_read = handle.read()
        handle.close()

        root = ET.fromstring(handle_read)

        pmcid = ""

        for link in root.iter('Link'):
            for id in link.iter('Id'):
                pmcid = id.text

        if pmcid != "":
            print(f"PMC ID for PMID {pmid}: {pmcid}")
        else:
            print(f"No PMC ID found for PMID {pmid}")

        return pmcid != "", pmcid if pmcid else None

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

        if pmid == "None":
            print(f"Full text not available in PMC for PMID {pmid}")
            return None

        print(f"Fetching full text for PMC ID {pmid}")
        content = ""
        retstart = 0

        while True:
            full_text_handle = Entrez.efetch(
                db="pmc", id=pmid, rettype="xml", retstart=retstart
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
            time.sleep(2)

        result = extract_text(content)
        return result

    except Exception as e:
        print(f"Error getting full text for PMID {pmid}: {str(e)}")
        return None


def sanitize_filename(name):
    """Sanitize the filename to remove/replace invalid characters."""
    name = str(name)
    name = re.sub(r'[\\/*?:"<>|]', "_", name)  # Windows-safe
    name = re.sub(r'\s+', '_', name)  # Replace spaces with underscores
    name = re.sub(r'[^\w\-_.]', '_', name)  # Only allow safe characters
    return name.strip('_')

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



def fetch_and_save(pmid, name, output_dir, errors_csv_file):
    """Worker function for downloading and saving full text."""
    print(f"Trying to fetch pmid {pmid}")
    safe_name = sanitize_filename(name)
    
    try:
        output = get_full_text(str(pmid))
        if output:
            filepath = os.path.join(output_dir, f"{safe_name}.txt")
            with open(filepath, "w", encoding="utf-8") as f:
                f.write(output)
            print(f"** fetching of reprint {pmid} succeeded")
        else:
            raise ValueError("No full text available")

    except requests.HTTPError as e:
        status = e.response.status_code
        print(f"** fetching of reprint {pmid} failed, HTTP {status}")
        with lock:
            with open(errors_csv_file, 'a', encoding='utf-8') as f:
                f.write(f"{pmid},{safe_name},HTTP {status}\n")

    except Exception as e:
        print(f"** fetching of reprint {pmid} failed with error: {str(e)}")
        with lock:
            with open(errors_csv_file, 'a', encoding='utf-8') as f:
                f.write(f"{pmid},{safe_name},{str(e)}\n")


def run_pubmed_downloader(output_dir, input_csv_file, errors_csv_file, max_workers=8):
    """Main execution function with multithreading."""
    if not os.path.exists(output_dir):
        print(f"Output directory {output_dir} did not exist. Created.")
        os.makedirs(output_dir)

    df = pd.read_csv(input_csv_file)
    pmids = df["PMCID"].tolist()
    names = df["Title"].astype(str).tolist()

    # Initialize error log file if it doesn't exist
    if not os.path.exists(errors_csv_file):
        with open(errors_csv_file, 'w', encoding='utf-8') as f:
            f.write("pmid,name,error\n")

    # Use ThreadPoolExecutor to run downloads in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(fetch_and_save, pmid, name, output_dir, errors_csv_file)
            for pmid, name in zip(pmids, names)
        ]
        # Optional: wait for all to complete
        concurrent.futures.wait(futures)

if __name__ == "__main__":
    args = parse_arguments()
    run_pubmed_downloader(args["out"], args["csv"], args["errors"])
