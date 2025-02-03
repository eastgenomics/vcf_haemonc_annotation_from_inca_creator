import pandas as pd
import argparse

# Load CSV into a DataFrame with argparse
parser = argparse.ArgumentParser()
parser.add_argument('-filename', '--filename')

args = parser.parse_args()

# Ensure 'date_last_evaluated' is parsed as datetime
df = pd.read_csv(args.filename, parse_dates=['date_last_evaluated'])

# Remove any spaces and substitute with _
df['oncogenicity_classification'] = df['oncogenicity_classification'].str.replace(' ', '_')
df['chromosome'] = df['chromosome'].str.replace(' ', '')

# Sort and flag the latest variant records
df_sorted = df.sort_values(
    by=['hgvsc', 'date_last_evaluated'],
    ascending=[True, False]
)

# Flag the latest record for each group
df_sorted['is_latest'] = df_sorted.gry(['hgvsc']).cumcount() == 0

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
print(merged_df.shape)

# Write the VCF file
with open('output.vcf', 'w') as vcf_file:
    # Write VCF header LC = latest classification
    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("##reference=GRCh38\n")
    vcf_file.write('##INFO=<ID=EG_LC,Number=1,Type=String,Description="Oncogenicity classification of the variant">\n')
    vcf_file.write('##INFO=<ID=EG_LCD,Number=1,Type=String,Description="Date of latest classification evaluation">\n')
    vcf_file.write('##INFO=<ID=EG_LCS,Number=1,Type=String,Description="Specimen ID of the latest classification">\n')
    vcf_file.write('##INFO=<ID=EG_OC,Number=.,Type=String,Description="Other Classifications with counts">\n')
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
            f"EG_OC={other_oncogenicity}"
        )
        vcf_file.write(
            f"{row['chromosome']}\t{row['start']}\t{row['hgvsc']}\t"
            f"{row['reference_allele']}\t{row['alternate_allele']}\t.\t.\t{info_field}\n"
        )

print("VCF file created as 'output.vcf'")