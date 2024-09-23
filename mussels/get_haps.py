import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
plt.rcParams.update({'font.size': 7.0})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
plt.rcParams['svg.fonttype'] = 'none'

sns.set_palette("muted")

# Initialize variables
#file_path = "./dbghaplo.mussel"
file_path = "./all-min5/haplotypes.fasta"
data = {}

# Regular expression to capture the header line fields
header_pattern = re.compile(r">Contig:(.*),Range:(.*),Haplotype:(.*),Abundance:(.*),Depth:(.*)")

# Parse the FASTA-like file
with open(file_path, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            match = header_pattern.match(line)
            if match:
                contig = match.group(1)
                range_ = match.group(2)
                
                # Use (contig, range) tuple as a key in the dictionary
                key = (contig, range_)
                
                # Increment the haplotype count for the given contig and range
                if key not in data:
                    data[key] = 0
                data[key] += 1

# Convert the dictionary to a pandas dataframe
df = pd.DataFrame([(k[0], k[1], v) for k, v in data.items()], columns=["Contig", "Range", "Haplotypes"])

# Print the dataframe
print(df)
fig,ax = plt.subplots(figsize=(2,2))
sns.set_palette("muted")
ax = sns.boxplot(df["Haplotypes"])
#ax = sns.swarmplot(y='Haplotypes', data=df, color="black")
ax.set_ylabel("Predicted # of haplotypes\nper region")
ax.set_title("dbghaplo haplotyping on\ncandidate regions")

#dashed line at 11
plt.axhline(y=11, color='black', linestyle='--', label="Published estimate (11 haplotypes)")

med = df["Haplotypes"].median()
plt.axhline(y=med, color='blue', linestyle='--', label=f"dbghaplo median ({int(med)} haplotypes)")
plt.axhline(y=2, color='red', linestyle='--', label="igda median (2 haplotypes)")
plt.axhline(y=1, color='orange', linestyle='--', label="rvhaplo (1 haplotypes)")

plt.legend(frameon=False, loc='upper right', fontsize=7)

#despine
sns.despine()

plt.savefig("../figures/mussel_boxplot.svg")
plt.show()
# Optionally, save the dataframe to a CSV file
# df.to_csv("haplotypes_summary.csv", index=False)

