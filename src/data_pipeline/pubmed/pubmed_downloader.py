#!/usr/bin/env python
"""
PubMed Central PDF Fetcher

This tool downloads PDF articles from various publishers using PubMed IDs.
It supports multiple publishers including NEJM, Science Direct, and PubMed Central.
"""

import os
import re
import sys
import argparse
import urllib.parse
from bs4 import BeautifulSoup
import requests
import pandas as pd

def get_main_url(url):
    """Extract the base URL from a full URL."""
    return "/".join(url.split("/")[:3])


def save_pdf_from_url(pdf_url, directory, name, headers):
    """Download and save a PDF file from a given URL."""
    try:
        response = requests.get(pdf_url, headers=headers, allow_redirects=True)
        response.raise_for_status()  # Raise exception for HTTP errors

        filepath = f'{directory}/{name}.pdf'
        with open(filepath, 'wb') as f:
            f.write(response.content)
        return True
    except Exception as e:
        print(f"Error saving PDF: {e}")
        return False


class PdfFinder:
    """Contains methods to locate PDF download links from various publishers."""

    @staticmethod
    def acs_publications(req, soup, headers):
        """Find PDF links on ACS Publications websites."""
        possibleLinks = [x for x in soup.find_all('a') if type(x.get('title')) == str and (
            'high-res pdf' in x.get('title').lower() or 'low-res pdf' in x.get('title').lower())]

        if possibleLinks:
            print("** fetching reprint using the 'acsPublications' finder...")
            return get_main_url(req.url) + possibleLinks[0].get('href')
        return None

    @staticmethod
    def direct_pdf_link(req, soup, headers):
        """Check if URL directly points to a PDF."""
        if req.content[-4:] == b'.pdf':
            print("** fetching reprint using the 'direct pdf link' finder...")
            return req.url
        return None

    @staticmethod
    def future_medicine(req, soup, headers):
        """Find PDF links on Future Medicine websites."""
        possibleLinks = soup.find_all(
            'a', attrs={'href': re.compile("/doi/pdf")})
        if possibleLinks:
            print("** fetching reprint using the 'future medicine' finder...")
            return get_main_url(req.url) + possibleLinks[0].get('href')
        return None

    @staticmethod
    def generic_citation_labelled(req, soup, headers):
        """Find PDF links using citation_pdf_url meta tag."""
        possibleLinks = soup.find_all(
            'meta', attrs={'name': 'citation_pdf_url'})
        if possibleLinks:
            print("** fetching reprint using the 'generic citation labelled' finder...")
            return possibleLinks[0].get('content')
        return None

    @staticmethod
    def nejm(req, soup, headers):
        """Find PDF links on New England Journal of Medicine website."""
        possibleLinks = [x for x in soup.find_all('a') if type(x.get(
            'data-download-type')) == str and (x.get('data-download-type').lower() == 'article pdf')]

        if possibleLinks:
            print("** fetching reprint using the 'NEJM' finder...")
            return get_main_url(req.url) + possibleLinks[0].get('href')
        return None

    @staticmethod
    def pubmed_central_v1(req, soup, headers):
        """Find PDF links on PubMed Central (version 1)."""
        possibleLinks = soup.find_all('a', re.compile('pdf'))
        # Filter out epdf links (for Wiley compatibility)
        possibleLinks = [
            x for x in possibleLinks if 'title' in x.attrs and 'epdf' not in x.get('title').lower()]

        if possibleLinks:
            print("** fetching reprint using the 'pubmed central v1' finder...")
            return get_main_url(req.url) + possibleLinks[0].get('href')
        return None

    @staticmethod
    def pubmed_central_v2(req, soup, headers):
        """Find PDF links on PubMed Central (version 2)."""
        possibleLinks = soup.find_all(
            'a', attrs={'href': re.compile('/pmc/articles')})

        if possibleLinks:
            print("** fetching reprint using the 'pubmed central v2' finder...")
            return f"https://www.ncbi.nlm.nih.gov{possibleLinks[0].get('href')}"
        return None

    @staticmethod
    def science_direct(req, soup, headers):
        """Find PDF links on Science Direct websites."""
        try:
            newUri = urllib.parse.unquote(
                soup.find_all('input')[0].get('value'))
            req = requests.get(newUri, allow_redirects=True, headers=headers)
            soup = BeautifulSoup(req.content, 'lxml')

            possibleLinks = soup.find_all(
                'meta', attrs={'name': 'citation_pdf_url'})

            if possibleLinks:
                print("** fetching reprint using the 'science_direct' finder...")
                req = requests.get(possibleLinks[0].get(
                    'content'), headers=headers)
                soup = BeautifulSoup(req.content, 'lxml')
                return soup.find_all('a')[0].get('href')
        except (IndexError, AttributeError):
            pass
        return None

    @staticmethod
    def uchicago_press(req, soup, headers):
        """Find PDF links on University of Chicago Press websites."""
        possibleLinks = [x for x in soup.find_all('a') if type(x.get(
            'href')) == str and 'pdf' in x.get('href') and '.edu/doi/' in x.get('href')]

        if possibleLinks:
            print("** fetching reprint using the 'uchicagoPress' finder...")
            return get_main_url(req.url) + possibleLinks[0].get('href')
        return None


def fetch(pmid, finders_list, name, headers, error_file):
    """
    Attempt to fetch a PDF for the given PubMed ID using various finder methods.
    
    Args:
        pmid (str): PubMed ID to fetch
        finders_list (list): List of finder method names to try
        name (str): Name to save the PDF as
        headers (dict): HTTP headers for requests
        error_file: File to write errors to
        
    Returns:
        bool: Success status
    """
    # Get the redirect URL from PubMed
    uri = f"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id={pmid}&retmode=ref&cmd=prlinks"

    try:
        req = requests.get(uri, headers=headers)
        req.raise_for_status()

        # Skip Ovid links as they're not supported
        if 'ovid' in req.url:
            print(
                f" ** Reprint {pmid} cannot be fetched as ovid is not supported by the requests package.")
            error_file.write(f"{pmid}\t{name}\n")
            return True

        soup = BeautifulSoup(req.content, 'lxml')

        # Try each finder method until one succeeds
        for finder_name in finders_list:
            print(f"Trying {finder_name}")
            finder_method = getattr(PdfFinder, finder_name)
            pdf_url = finder_method(req, soup, headers)

            if pdf_url is not None:
                if save_pdf_from_url(pdf_url, args['out'], name, headers):
                    print(f"** fetching of reprint {pmid} succeeded")
                    return True

        # If we get here, all finders failed
        print(
            f"** Reprint {pmid} could not be fetched with the current finders.")
        error_file.write(f"{pmid}\t{name}\n")
        return False

    except Exception as e:
        print(f"** Error in fetch function: {e}")
        error_file.write(f"{pmid}\t{name}\n")
        return False


def parse_arguments():
    """Parse and validate command line arguments."""
    parser = argparse.ArgumentParser()
    parser._optionals.title = "Flag Arguments"
    parser.add_argument(
        '-pmids', help="Comma separated list of pmids to fetch. Must include -pmids or -pmf.", default='%#$')
    parser.add_argument(
        '-csv', help="File with pmids to fetch inside, one pmid per line. Optionally, the file can be a tsv with a second column of names to save each pmid's article with (without '.pdf' at the end). Must include -pmids or -pmf", default='%#$')
    parser.add_argument(
        '-out', help="Output directory for fetched articles. Default: fetched_pdfs", default="fetched_pdfs")
    parser.add_argument(
        '-errors', help="Output file path for pmids which failed to fetch. Default: unfetched_pmids.tsv", default="unfetched_pmids.tsv")
    parser.add_argument(
        '-maxRetries', help="Change max number of retries per article on an error 104. Default: 3", default=3, type=int)

    args = vars(parser.parse_args())

    # Check for required arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args['pmids'] == '%#$' and args['csv'] == '%#$':
        print("Error: Either -pmids or -pmf must be used. Exiting.")
        sys.exit(1)

    if args['pmids'] != '%#$' and args['csv'] != '%#$':
        print("Error: -pmids and -pmf cannot be used together. Ignoring -pmf argument")
        args['pmf'] = '%#$'

    return args


def run_pubmed_downloader(output_dir, input_csv_file, errors_csv_file, max_retries=3):
    """Main execution function."""
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        print(
            f"Output directory of {output_dir} did not exist. Created the directory.")
        os.mkdir(output_dir)

    # Set up headers for HTTP requests
    headers = requests.utils.default_headers()
    headers['User-Agent'] = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/56.0.2924.87 Safari/537.36'

    # Set up list of finder methods (converted to snake_case)
    finders = [
        'generic_citation_labelled',
        'pubmed_central_v2',
        'acs_publications',
        'uchicago_press',
        'nejm',
        'future_medicine',
        'science_direct',
        'direct_pdf_link',
    ]

    df = pd.read_csv(input_csv_file)
    pmids = df['PMID'].tolist()
    names = df['Title'].tolist()
    names = [str(name) for name in names]  # Ensure all names are strings

    names = [name.replace(" ", "_").replace(" ", "_").replace(".", "_") for name in names]

    # Process each PMID
    errors = []
    
    for pmid, name in zip(pmids, names):
        print(f"Trying to fetch pmid {pmid}")
        
        retries = 0
        while retries < max_retries:
            try:
                fetch(pmid, finders, name, headers, None)
                break
            except requests.ConnectionError as e:
                if '104' in str(e) or 'BadStatusLine' in str(e):
                    retries += 1
                    if retries < max_retries:
                        print(f"** fetching of reprint {pmid} failed from error {e}, retrying")
                    else:
                        print(f"** fetching of reprint {pmid} failed from error {e}")
                        errors.append({'pmid': pmid, 'name': name})
                else:
                    print(f"** fetching of reprint {pmid} failed from error {e}")
                    errors.append({'pmid': pmid, 'name': name})
                    break
            except Exception as e:
                print(f"** fetching of reprint {pmid} failed from error {e}")
                errors.append({'pmid': pmid, 'name': name})
                break
    
    # Write errors to CSV
    if errors:
        pd.DataFrame(errors).to_csv(errors_csv_file, index=False, header=False)


if __name__ == "__main__":
    args = parse_arguments()
    run_pubmed_downloader(args['out'], args['csv'], args['errors'], args['maxRetries'])
