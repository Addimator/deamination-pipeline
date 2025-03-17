import pandas as pd
import altair as alt

# file_path = "/home/adrian/Documents/Promotion/deamination_stuff/deamination-pipeline/results/debug/number_bases.csv"
file_path = snakemake.input[0]
df = pd.read_csv(file_path)

# Daten umstrukturieren für Altair
df_melted = df.melt(
    id_vars=["direction", "methylation_status", "orig_base", "abs_number"],
    value_vars=["A", "C", "G", "T"],
    var_name="Nucleotide",
    value_name="Count",
)

df_melted["ratio"] = (df_melted["Count"] / df_melted["abs_number"]).round(2)

# Kategorien für die Gruppierung definieren
categories = [
    (False, "forward"),
    (True, "forward"),
    (False, "reverse"),
    (True, "reverse"),
]

charts = []
for methylation, direction in categories:
    subset = df_melted[
        (df_melted["methylation_status"] == methylation)
        & (df_melted["direction"] == direction)
    ]
    print(subset)

    chart = (
        alt.Chart(
            subset,
        )
        .mark_bar()
        .encode(
            x=alt.X("Nucleotide:N", title="Nucleotide"),
            y=alt.Y("Count:Q", title="Count"),
            color="Nucleotide:N",
            column=alt.Column("orig_base:N", title="CpG Position"),
            row=alt.Row("methylation_status:N", title="Methylation"),
            tooltip=["Count:N", "ratio:N"],
        )
        .properties(
            width=100,
            height=150,
        )
    )

    charts.append(chart)

# Raster als 2×4 Anordnung (zwei Zeilen mit vier Plots)
final_chart = alt.vconcat(
    alt.hconcat(*charts[:2]).properties(
        title=alt.TitleParams(
            text="Forward Strand",
        )
    ),
    alt.hconcat(*charts[2:]).properties(title=alt.TitleParams(text="Reverse Strand")),
)

# Speichern als HTML
output_path = (
    # "/home/adrian/Documents/Promotion/deamination_stuff/deamination-pipeline/test.html"
    snakemake.output[0]
)
final_chart.save(output_path, scale_factor=2.0)
