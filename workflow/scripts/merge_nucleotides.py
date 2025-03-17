import csv
from collections import defaultdict
import pandas as pd
import altair as alt
from collections import defaultdict


def process_csv(file_path):
    # CSV einlesen
    df = pd.read_csv(file_path)

    # Sicherstellen, dass numerische Spalten als Integer behandelt werden
    df["orig_base"] = df["orig_base"]
    for base in ["A", "C", "G", "T", "N"]:
        df[base] = df[base].astype(int)

    # Gruppieren nach 'direction', 'methylation_status' und 'orig_base'
    grouped_df = df.groupby(["methylation_status", "direction", "orig_base"])[
        ["A", "C", "G", "T", "N"]
    ].sum()

    grouped_df["abs_number"] = grouped_df[["A", "C", "G", "T", "N"]].sum(axis=1)

    return grouped_df.reset_index()


# Beispielaufruf
# file_path = "/home/adrian/Documents/Promotion/deamination_stuff/deamination-pipeline/results/debug/number_bases_concat.csv"
file_path = snakemake.input[0]

result_df = process_csv(file_path)
print(result_df)

# out_path = "/home/adrian/Documents/Promotion/deamination_stuff/deamination-pipeline/results/debug/number_bases.csv"
out_path = snakemake.output[0]

result_df.to_csv(out_path, index=False)
