import altair as alt
import pandas as pd
import numpy as np

import heapq


meth_pos = []
unmeth_pos = []
forward_dict = {}
reverse_dict = {}

# Create dictionaries for Bedgraph and VCF
with open(snakemake.input["bedGraph"], 'r') as ref_file, open(snakemake.input["ref_bases"][0], 'r') as sd_file:

    for line in sd_file:
        if not line.startswith('#'):
            parts = line.strip().split('\t')
            chrom, pos, dir, a, c, g, t, n = parts[0], int(
                parts[1]), parts[2], *map(int, parts[3:])
            if dir == "f_0":
                forward_dict[(chrom, pos)] = t / (a + c + g + n)
            elif dir == "r_1":
                reverse_dict[(chrom, pos)] = a / (c + g + t + n)

    for line in ref_file:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value = parts[0].replace(
                "chr", ""), (int(parts[1]) + int(parts[2])) // 2, float(parts[3])
            if methylation_value > 0:
                meth_pos.append(chrom, position)
            else:
                unmeth_pos.append(chrom, position)


meth_forward = {pos: forward_dict[pos]
                for pos in meth_pos if pos in forward_dict}
meth_reverse = {pos: reverse_dict[pos]
                for pos in meth_pos if pos in reverse_dict}
unmeth_forward = {pos: forward_dict[pos]
                  for pos in unmeth_pos if pos in forward_dict}
unmeth_reverse = {pos: reverse_dict[pos]
                  for pos in unmeth_pos if pos in reverse_dict}

sorted_meth_forward = dict(
    sorted(meth_forward.items(), key=lambda item: item[0][1]))
sorted_meth_reverse = dict(
    sorted(meth_reverse.items(), key=lambda item: item[0][1]))
sorted_unmeth_forward = dict(
    sorted(unmeth_forward.items(), key=lambda item: item[0][1]))
sorted_unmeth_reverse = dict(
    sorted(unmeth_reverse.items(), key=lambda item: item[0][1]))


# Plot methylated positions
numbers_list = list(
    range(1, max(len(sorted_meth_forward), len(sorted_meth_reverse)) + 1))

data_forward = pd.DataFrame({
    'Methylation': [sorted_meth_forward.get(pos, None) for pos in numbers_list],
    'Positions': numbers_list,
    'Line': ['Forward'] * len(numbers_list)
})

# Daten für die zweite Linie (sorted_meth_reverse)
data_reverse = pd.DataFrame({
    'Methylation': [sorted_meth_reverse.get(pos, None) for pos in numbers_list],
    'Positions': numbers_list,
    'Line': ['Reverse'] * len(numbers_list)
})

# Kombiniere die Daten für beide Linien
data_combined = pd.concat([data_forward, data_reverse])

# Erstelle den Plot
scatter = alt.Chart(data_combined).mark_circle(opacity=0.5).encode(
    x='Methylation',
    y='Positions',
    color='Line:N'
)

final_chart = scatter.properties(
    width=400,
    height=400,
    title='Methylated positions'
)

final_chart.save(snakemake.output["meth"], scale_factor=2.0)


# Plot unmethylated positions
numbers_list = list(
    range(1, max(len(sorted_unmeth_forward), len(sorted_unmeth_reverse)) + 1))

data_forward = pd.DataFrame({
    'Unmethylation': [sorted_unmeth_forward.get(pos, None) for pos in numbers_list],
    'Positions': numbers_list,
    'Line': ['Forward'] * len(numbers_list)
})

# Daten für die zweite Linie (sorted_unmeth_reverse)
data_reverse = pd.DataFrame({
    'Unethylation': [sorted_unmeth_reverse.get(pos, None) for pos in numbers_list],
    'Positions': numbers_list,
    'Line': ['Reverse'] * len(numbers_list)
})

# Kombiniere die Daten für beide Linien
data_combined = pd.concat([data_forward, data_reverse])

# Erstelle den Plot
scatter = alt.Chart(data_combined).mark_circle(opacity=0.5).encode(
    x='Unmethylation',
    y='Positions',
    color='Line:N'
)

final_chart = scatter.properties(
    width=400,
    height=400,
    title='Unmethylated positions'
)

final_chart.save(snakemake.output["unmeth"], scale_factor=2.0)
