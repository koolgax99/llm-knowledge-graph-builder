import pandas as pd
import requests
import time

def convert_pmids_to_pmcids_with_full_row(input_csv_path, output_csv_path, tool_name, email):
    """
    Maps PMIDs to PMCIDs and retains the full input row in the output.

    Parameters:
        input_csv_path (str): CSV with a column named "PMID".
        output_csv_path (str): Destination CSV with all original columns + PMCID.
        tool_name (str): Application name for NCBI API.
        email (str): Contact email address for NCBI API.
    """
    # Load original data
    df = pd.read_csv(input_csv_path, dtype={"PMID": str})
    pmid_list = df["PMID"].tolist()

    API_URL = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
    results = []

    # Fetch PMCID for batches of PMIDs
    for i in range(0, len(pmid_list), 200):
        batch = pmid_list[i:i + 200]
        params = {
            "tool": tool_name,
            "email": email,
            "ids": ",".join(batch),
            "idtype": "pmid",
            "format": "json"
        }
        response = requests.get(API_URL, params=params)
        response.raise_for_status()
        data = response.json()
        for rec in data.get("records", []):
            results.append({
                "PMID": str(rec.get("pmid")),
                "PMCID": rec.get("pmcid", "None")
            })
        time.sleep(2)  # NCBI throttle

    # Create mapping DataFrame and merge with original
    mapping_df = pd.DataFrame(results)
    merged_df = df.merge(mapping_df, on="PMID", how="left")

    # Save full result with PMCID appended
    merged_df.to_csv(output_csv_path, index=False)
    print(f"✅ Full row output saved to '{output_csv_path}'")

# Example usage:
# convert_pmids_to_pmcids_with_full_row(
#     input_csv_path="pubmed_plasticity_pubmed_results.csv",
#     output_csv_path="pmid_to_pmcid_full_output.csv",
#     tool_name="nihar_pmcid_mapper",
#     email="sanda.n@northeastern.edu"
# )
