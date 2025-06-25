import pandas as pd
import requests
import os

# Supplementary table 2 from (Mathys 2019, Nature) --> https://www.nature.com/articles/s41586-019-1195-2
URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1195-2/MediaObjects/41586_2019_1195_MOESM4_ESM.xlsx"
FILENAME = "data/mathys2019_supplementary_table2.xlsx"
CELLTYPES = ["Ex", "In", "Ast", "Oli", "Opc", "Mic"]  # Sheets in the Excel file

if not os.path.exists(FILENAME):  # Download the file if we don't have it yet.
    print(f"Downloading {URL} to {FILENAME}...")
    response = requests.get(URL)
    response.raise_for_status()
    with open(FILENAME, "wb") as f:
        f.write(response.content)

ct_dfs = {}
degs = set()
for sheet in CELLTYPES:
    print(f"\nProcessing sheet: {sheet}")
    # skip the first row
    # only use cols A-I
    # use the first column as index
    # use the columns as they are
    df = pd.read_excel(FILENAME, sheet_name=sheet, header=1, usecols="A:I", index_col=0)
    
    DEG_COL = "DEGs.Ind.Mix.models"
    df = df[df[DEG_COL]]  # Filter out rows where the DEG column is NaN

    print(df.head())
    print(f"Found {len(df)} DEGs in {sheet} cell type.")

    ct_dfs[sheet] = df
    
    degs.update(df.index)  # Collect all genes across all cell types

print(f"\nFound {len(degs)} unique DEGs across all cell types.")

# Combine everything into a single DataFrame and 
combined_df = pd.concat(ct_dfs, axis=1)
combined_df.to_csv("data/mathys2019_DEGs.csv")
with open("data/mathys2019_DEGs_genes.txt", "w") as f:
    f.write("\n".join(list(degs)))
