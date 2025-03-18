import urllib.request
import time
import os


bedGraphs_params = snakemake.params["bedGraphs"]
file_path = snakemake.output[0]

file_no_gz, _ = os.path.splitext(
    file_path
)  # "resources/HG002/GSM5649436_TruSeq_HG002_LAB01_REP02.bedGraph"

# Dann ".bedGraph" entfernen
bedGraph_name, _ = os.path.splitext(
    file_no_gz
)  # "resources/HG002/GSM5649436_TruSeq_HG002_LAB01_REP02"

# Nur den Dateinamen ohne Ordner
bedGraph = os.path.basename(bedGraph_name)


print(bedGraphs_params)
print("test")
print(bedGraph)
print("test")


def download_with_retries(url, output_path, retries=5, delay=5):
    for attempt in range(retries):
        try:
            print(f"Attempt {attempt+1}/{retries}: Downloading {url}")
            urllib.request.urlretrieve(url, output_path)
            print("Download successful!")
            return  # Exit the function if successful
        except urllib.error.ContentTooShortError as e:
            print(
                f"Download failed (ContentTooShortError), retrying in {delay} seconds..."
            )
        except urllib.error.URLError as e:
            print(f"Download failed (URLError: {e}), retrying in {delay} seconds...")
        except Exception as e:
            print(f"Unexpected error: {e}, retrying in {delay} seconds...")

        time.sleep(delay)

    print(f"Failed to download {url} after {retries} attempts.")


# for bedGraph in bedGraphs:
accession_number = bedGraph.split("_")[0]
url = f"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5649nnn/{accession_number}/suppl/{bedGraph}.bedGraph.gz"
output_path = f"resources/HG002/{bedGraph}.bedGraph.gz"
download_with_retries(url, output_path)
