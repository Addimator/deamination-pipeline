import altair as alt
import pandas as pd
import numpy as np

# alt.data_transformers.disable_max_rows()
alt.data_transformers.enable("vegafusion")
meth_pos = []
unmeth_pos = []
forward_dict = {}
reverse_dict = {}
print(snakemake.input["bedGraph"])
print(snakemake.input["ref_bases"])
# Create dictionaries for Bedgraph and VCF
with open(snakemake.input["bedGraph"][0], "r") as ref_file, open(
    snakemake.input["ref_bases"], "r"
) as sd_file:

    for line in sd_file:
        if not line.startswith("#"):
            parts = line.strip().split("\t")
            chrom, pos, dir, a, c, g, t, n = (
                parts[0],
                int(parts[1]),
                parts[2],
                *map(int, parts[3:]),
            )
            if dir == "f_0":
                forward_dict[(chrom, pos)] = (
                    1 if (a + c + g + t + n) == 0 else t / (a + c + g + t + n)
                )

            elif dir == "r_1":
                reverse_dict[(chrom, pos)] = (
                    1 if (a + c + g + t + n) == 0 else a / (a + c + g + t + n)
                )

    for line in ref_file:
        if not line.startswith("track"):
            parts = line.strip().split("\t")
            chrom, position, methylation_value = (
                parts[0].replace("chr", ""),
                (int(parts[1]) + int(parts[2])) // 2,
                float(parts[3]),
            )
            if methylation_value > 0:
                meth_pos.append((chrom, position))
            else:
                unmeth_pos.append((chrom, position))


meth_forward = {pos: forward_dict[pos] for pos in meth_pos if pos in forward_dict}
meth_reverse = {pos: reverse_dict[pos] for pos in meth_pos if pos in reverse_dict}
unmeth_forward = {pos: forward_dict[pos] for pos in unmeth_pos if pos in forward_dict}
unmeth_reverse = {pos: reverse_dict[pos] for pos in unmeth_pos if pos in reverse_dict}

sorted_meth_forward = dict(sorted(meth_forward.items(), key=lambda item: item[0][1]))
sorted_meth_reverse = dict(sorted(meth_reverse.items(), key=lambda item: item[0][1]))
sorted_unmeth_forward = dict(
    sorted(unmeth_forward.items(), key=lambda item: item[0][1])
)
sorted_unmeth_reverse = dict(
    sorted(unmeth_reverse.items(), key=lambda item: item[0][1])
)


# Plot methylated positions
# numbers_list = list(
#     range(1, max(len(sorted_meth_forward), len(sorted_meth_reverse)) + 1))
# Daten für die erste Linie (sorted_meth_forward)
data_forward = pd.DataFrame(
    {
        "Methylation": list(sorted_meth_forward.values()),
        "Positions": list(sorted_meth_forward.keys()),
        "Line": ["Forward"] * len(sorted_meth_forward),
    }
)

# Daten für die zweite Linie (sorted_meth_reverse)
data_reverse = pd.DataFrame(
    {
        "Methylation": list(sorted_meth_reverse.values()),
        "Positions": list(sorted_meth_reverse.keys()),
        "Line": ["Reverse"] * len(sorted_meth_reverse),
    }
)

# Kombiniere die Daten für beide Linien
data_combined = pd.concat([data_forward, data_reverse])


# Zähle die Anzahl der Zeilen mit Methylation = 0.0
count_zero_methylation = data_combined[data_combined["Methylation"] == 0.0].shape[0]
count_one_methylation = data_combined[data_combined["Methylation"] != 0.0].shape[0]
print(f"Anzahl der Zeilen mit Methylation = 0.0: {count_zero_methylation}")
print(f"Anzahl der Zeilen mit Methylation = 1.0: {count_one_methylation}")

# Entferne alle Zeilen mit Methylation = 0.0
df_filtered = data_combined[data_combined["Methylation"] != 0.0]


# Erstelle den Plot
scatter = (
    alt.Chart(df_filtered)
    .mark_circle(opacity=0.5)
    .encode(y="Methylation", x="Positions", color="Line:N")
)

final_chart = scatter.properties(width=400, height=400, title="Methylated positions")

final_chart.save(snakemake.output["meth"], scale_factor=2.0)


# Plot unmethylated positions

data_forward = pd.DataFrame(
    {
        "Unmethylation": sorted_unmeth_forward.values(),
        "Positions": sorted_unmeth_forward.keys(),
        "Line": ["Forward"] * len(sorted_unmeth_forward),
    }
)

# Daten für die zweite Linie (sorted_unmeth_reverse)
data_reverse = pd.DataFrame(
    {
        "Unmethylation": sorted_unmeth_reverse.values(),
        "Positions": sorted_unmeth_reverse.keys(),
        "Line": ["Reverse"] * len(sorted_unmeth_reverse),
    }
)

# Kombiniere die Daten für beide Linien
data_combined = pd.concat([data_forward.head(2000), data_reverse.head(2000)])


count_zero_methylation = data_combined[data_combined["Unmethylation"] == 0.0].shape[0]
print(f"Anzahl der Zeilen mit Unmethylation = 0.0: {count_zero_methylation}")

# Entferne alle Zeilen mit Methylation = 0.0
df_filtered = data_combined[data_combined["Unmethylation"] != 0.0]


# Erstelle den Plot
scatter = (
    alt.Chart(df_filtered)
    .mark_circle(opacity=0.5)
    .encode(y="Unmethylation", x="Positions", color="Line:N")
)

final_chart = scatter.properties(width=400, height=400, title="Unmethylated positions")

final_chart.save(snakemake.output["unmeth"], scale_factor=2.0)
