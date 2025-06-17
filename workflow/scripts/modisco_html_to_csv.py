import pandas as pd
from bs4 import BeautifulSoup
import os

def extract_table_from_html(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        soup = BeautifulSoup(f, "html.parser")
    
    table = soup.find("table", class_="dataframe")
    if table is None:
        return None
    
    headers = [th.text for th in table.find("thead").find_all("th")]
    rows = []
    for tr in table.find("tbody").find_all("tr"):
        cells = [td.text.strip() for td in tr.find_all("td")]
        rows.append(cells)
    
    return pd.DataFrame(rows, columns=headers)

df = extract_table_from_html(snakemake.input['html'])

## add columns
df['cluster'] = snakemake.params['cluster']
df['fold'] = snakemake.params['fold']

## left 
cols_left = ['cluster', 'fold']
idx_cols = [c for c in df.columns if c in cols_left]
idx_cols += [c for c in df.columns if c not in cols_left and 'logo' not in c and 'cwm' not in c]
df = df[idx_cols]

# Save to CSV
df.to_csv(snakemake.output['csv'], sep=",", index=False)
