import pandas as pd
import argparse
import subprocess

# Load CSV name from command line with argparse
parser = argparse.ArgumentParser()
parser.add_argument('-filename', '--filename')
args = parser.parse_args()

# Load CSV from Inca into a DataFrame 
# Ensure 'date_last_evaluated' is parsed as datetime and reduce the risk of type inference errors
df = pd.read_csv(args.filename, parse_dates=['date_last_evaluated'], low_memory=False) 

# Remove spaces and clean fields
df['oncogenicity_classification'] = df['oncogenicity_classification'].str.replace(' ', '_')
df['chromosome'] = df['chromosome'].str.replace(' ', '')

# Sort and flag the latest variant records
df_sorted = df.sort_values(
    by=['hgvsc', 'date_last_evaluated'],
    ascending=[True, False]
)

# Flag the latest record for each group
df_sorted['is_latest'] = df_sorted.groupby(['hgvsc']).cumcount() == 0

# Separate the latest records
latest_df = df_sorted[df_sorted['is_latest']].drop(columns='is_latest')

# Separate the other records
other_df = df_sorted[~df_sorted['is_latest']].drop(columns='is_latest')


# Aggregate other oncogenicity classifications with counts
def format_oncogenicity_counts(x):
    counts = x.value_counts()
    formatted_counts = [f"{cls}({count})" for cls, count in counts.items()]
    return "|".join(formatted_counts)


grouped_other_oncogenicity = other_df.groupby(['hgvsc']).agg({
    'oncogenicity_classification': format_oncogenicity_counts
}).rename(columns={'oncogenicity_classification': 'Other_Classifications'}).reset_index()


# Merge the latest records with aggregated other classifications
merged_df = pd.merge(
    latest_df,
    grouped_other_oncogenicity,
    on=['hgvsc'],
    how='left'
)


# Define contigs
contigs = """##contig=<ID=1,length=248956422,assembly=hg38>
##contig=<ID=2,length=242193529,assembly=hg38>
##contig=<ID=3,length=198295559,assembly=hg38>
##contig=<ID=4,length=190214555,assembly=hg38>
##contig=<ID=5,length=181538259,assembly=hg38>
##contig=<ID=6,length=170805979,assembly=hg38>
##contig=<ID=7,length=159345973,assembly=hg38>
##contig=<ID=8,length=145138636,assembly=hg38>
##contig=<ID=9,length=138394717,assembly=hg38>
##contig=<ID=10,length=133797422,assembly=hg38>
##contig=<ID=11,length=135086622,assembly=hg38>
##contig=<ID=12,length=133275309,assembly=hg38>
##contig=<ID=13,length=114364328,assembly=hg38>
##contig=<ID=14,length=107043718,assembly=hg38>
##contig=<ID=15,length=101991189,assembly=hg38>
##contig=<ID=16,length=90338345,assembly=hg38>
##contig=<ID=17,length=83257441,assembly=hg38>
##contig=<ID=18,length=80373285,assembly=hg38>
##contig=<ID=19,length=58617616,assembly=hg38>
##contig=<ID=20,length=64444167,assembly=hg38>
##contig=<ID=21,length=46709983,assembly=hg38>
##contig=<ID=22,length=50818468,assembly=hg38>
##contig=<ID=X,length=156040895,assembly=hg38>
##contig=<ID=Y,length=57227415,assembly=hg38>
"""

# Write the VCF file
with open('output.vcf', 'w') as vcf_file:
    # Write VCF header
    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("##reference=GRCh38\n")
    vcf_file.write(contigs)  # Write contigs
    vcf_file.write('##INFO=<ID=EG_LC,Number=1,Type=String,Description="Oncogenicity classification of the latest variant">\n')
    vcf_file.write('##INFO=<ID=EG_LCD,Number=1,Type=String,Description="Date of latest classification evaluation">\n')
    vcf_file.write('##INFO=<ID=EG_LCS,Number=1,Type=String,Description="Specimen ID of the latest classification">\n')
    vcf_file.write('##INFO=<ID=EG_OC,Number=.,Type=String,Description="Other Classifications with counts">\n')
    vcf_file.write('##INFO=<ID=EG_HGVS,Number=1,Type=String,Description="HGVS notation of the variant">\n')
    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # Write VCF data
    for _, row in merged_df.iterrows():
        # Handle possible missing Other_Classifications
        other_oncogenicity = (
            row['Other_Classifications']
            if pd.notnull(row['Other_Classifications']) else '.'
        )
        info_field = (
            f"EG_LC={row['oncogenicity_classification']};"
            f"EG_LCD={row['date_last_evaluated'].strftime('%Y-%m-%d')};"
            f"EG_LCS={row['specimen_id']};"
            f"EG_OC={other_oncogenicity};"
            f"EG_HGVS={row['hgvsc']}"
        )
        vcf_file.write(
            f"{row['chromosome']}\t{row['start']}\t{row['local_id']}\t"
            f"{row['reference_allele']}\t{row['alternate_allele']}\t.\t.\t{info_field}\n"
        )

## Sort VCF
with open("sorted.vcf", "w") as sorted_vcf:
    subprocess.run(["vcf-sort", "output.vcf"], stdout=sorted_vcf, check=True)

print("VCF file created as 'sorted.vcf'")
