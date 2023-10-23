import pickle
import numpy as np
import altair as alt
import pandas as pd
import random

ref_bases_file = snakemake.input.ref_bases
bedGraph_file = snakemake.input.bedGraph
# bedGraph_files = [snakemake.input[i] for i in range(len(snakemake.input["bedgraphs"]))]


dackel_dict = {}
pos_to_bases_dict = {}
true_dict = {}

# Create dictionaries for Bedgraph and VCF
with open(ref_bases_file, 'r') as ref_bases, open(bedGraph_file, 'r') as bedGraph:

    for line in bedGraph:
        if not line.startswith("track"):
            parts = line.strip().split('\t')
            chrom, position, methylation_value, coverage = parts[0], (int(parts[1]) + int(parts[2])) // 2, int(parts[3]), int(parts[4]) + int(parts[5])
            if coverage > 10:
                dackel_dict[(chrom, position)] = methylation_value


    for line in ref_bases:
        if not line.startswith('#'):
            parts = line.strip().split('\t')
            chrom, pos = parts[0], int(parts[1])
            dir, num_a, num_t, cov = parts[2], float(parts[3]), float(parts[6]), int(parts[3]) + int(parts[4]) + int(parts[5]) + int(parts[6]) + int(parts[7])
            if dir == "f" and cov > 10:        
                pos_to_bases_dict[(chrom, pos)] = num_t / cov
            if dir == "r" and cov > 10:        
                pos_to_bases_dict[(chrom, pos)] = num_a / cov
            



bedgraph_positions = [key for key in dackel_dict if key in pos_to_bases_dict]
bedgraph_meth_values = [dackel_dict[key] for key in bedgraph_positions]

pos_to_bases_positions = [key for key in pos_to_bases_dict if key in dackel_dict]
pos_to_bases_values = [pos_to_bases_dict[key] * 100 for key in pos_to_bases_positions]


missing_positions1 = [key for key in dackel_dict if key not in pos_to_bases_dict]
missing_positions2 = [key for key in pos_to_bases_dict if key not in dackel_dict]



# with open(snakemake.output["test"], "w") as datei:
#     for i, el in enumerate(true_meth_values):
#         if abs(el - pos_to_bases_values[i]) > 20: 
#             datei.write(str(true_positions[i]) + "\t" + str(el) + "\n")
#             datei.write("Bedgraph_value:" + str(bedgraph_meth_values[i]) + "\t Callsvcf:  " + str(pos_to_bases_values[i]) + "\n\n")

############################################################################################################################################################


line = alt.Chart(pd.DataFrame({'x': [1, 100], 'y': [1, 100]})).mark_line(color='red').encode(
    x='x:Q',
    y='y:Q'
)

# Plot TrueMeth vs Varlociraptor
deviation = sum(abs(x - y) for x, y in zip(bedgraph_meth_values, pos_to_bases_values))

data = pd.DataFrame({
    'Methylationrate Dackel': bedgraph_meth_values,
    'C to T rate': pos_to_bases_values
})

scatter = alt.Chart(data).mark_circle(opacity=0.5).encode(
    x='Methylationrate Dackel:Q',
    y='C to T rate:Q'
)

final_chart = (scatter + line).properties(
    width=400,
    height=400, 
    title=f'C to T rates (Deviation: {deviation})'
)
final_chart.save(snakemake.output.rates, scale_factor=2.0) 
