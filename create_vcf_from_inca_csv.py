import pandas as pd
import argparse
import datetime

# Function to format total classifications with counts
def format_total_classifications(classifications):
    """
    Counts all classifications, including the latest classification.
    Returns classifications in the format: classification(count)|classification(count)
    """
    counts = classifications.value_counts()
    formatted_counts = [f"{classification}({count})" for classification, count in counts.items()]
    return "|".join(formatted_counts)

# Function to aggregate HGVS annotations
def aggregate_hgvs(hgvs_series):
    """
    Aggregates all unique HGVS
    Returns them joined by |
    """
    unique_hgvs = hgvs_series.dropna().unique()
    return "|".join(unique_hgvs)

# Load CSV name from command line with argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', required=True, help="Input CSV file containing variant data")
args = parser.parse_args()

# Load CSV into a DataFrame with error handling
try:
    df = pd.read_csv(
        args.filename,
        parse_dates=['date_last_evaluated'],
        low_memory=False
    )
    required_columns = [
        'date_last_evaluated',
        'oncogenicity_classification',
        'chromosome',
        'hgvsc'
    ]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")
except Exception as e:
    print(f"Error loading CSV file: {e}")
    exit(1)

# Clean and preprocess fields
df['oncogenicity_classification'] = df['oncogenicity_classification'].str.replace(' ', '_')
df['chromosome'] = df['chromosome'].str.replace(' ', '')

# Ensure 'date_last_evaluated' is not missing
df = df.dropna(subset=['date_last_evaluated'])

# Group by unique variant identifiers
grouped = df.groupby(['chromosome', 'start', 'reference_allele', 'alternate_allele'])

# Function to get the latest entry based on 'date_last_evaluated'
def get_latest_entry(sub_df):
    """
    Function that returns the latest evaluated date
    by using idx.max()
    """
    latest_idx = sub_df['date_last_evaluated'].idxmax()
    latest_entry = sub_df.loc[latest_idx]
    return latest_entry

# Aggregate data for each unique variant
aggregated_data = []
for _, group in grouped:
    latest_entry = get_latest_entry(group)
    latest_classification = latest_entry['oncogenicity_classification']
    latest_date = latest_entry['date_last_evaluated']
    latest_sample_id = latest_entry['specimen_id']
    hgvs = aggregate_hgvs(group['hgvsc'])
    total_classifications = format_total_classifications(group['oncogenicity_classification'])  # Now includes all classifications
    
    aggregated_data.append({
        'chromosome': latest_entry['chromosome'],
        'start': latest_entry['start'],
        'reference_allele': latest_entry['reference_allele'],
        'alternate_allele': latest_entry['alternate_allele'],
        'latest_classification': latest_classification,
        'latest_date': latest_date,
        'latest_sample_id': latest_sample_id,
        'total_classifications': total_classifications,  # Replacing OC with TC
        'hgvs': hgvs
    })

aggregated_df = pd.DataFrame(aggregated_data)

# Generate a timestamp for the output filename
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
output_filename = f"haemonc_annotation_VCF_{timestamp}.vcf"


# Sort the merged DataFrame by Chromosome and start position:
# Define chromosome order: numeric first, then X and Y
chromosome_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
# Convert 'chromosome' column to categorical type with specified order
aggregated_df['chromosome'] = pd.Categorical(aggregated_df['chromosome'], categories=chromosome_order, ordered=True)
# Sort first by chromosome order, then by start position
aggregated_df = aggregated_df.sort_values(by=['chromosome', 'start'])


# Write the VCF file
with open(output_filename, 'w') as vcf_file:
    # Write VCF header
    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("##source=EG\n")
    vcf_file.write("##reference=GRCh38\n")
    vcf_file.write('##INFO=<ID=LC,Number=1,Type=String,Description="Latest classification">\n')
    vcf_file.write('##INFO=<ID=LCD,Number=1,Type=String,Description="Latest classification date">\n')
    vcf_file.write('##INFO=<ID=LCS,Number=1,Type=String,Description="Latest classification sample ID">\n')
    vcf_file.write('##INFO=<ID=TC,Number=.,Type=String,Description="Total classifications with counts (including latest)">\n')
    vcf_file.write('##INFO=<ID=HGVS,Number=.,Type=String,Description="HGVS annotations">\n')
    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # Write VCF data
    for _, row in aggregated_df.iterrows():
        info_field = (
            f"LC={row['latest_classification']};"
            f"LCD={row['latest_date'].strftime('%Y-%m-%d')};"
            f"LCS={row['latest_sample_id']};"
            f"TC={row['total_classifications']};"  # Changed OC to TC
            f"HGVS={row['hgvs']}"
        )
        vcf_file.write(
            f"{row['chromosome']}\t{row['start']}\t.\t"
            f"{row['reference_allele']}\t{row['alternate_allele']}\t.\t.\t{info_field}\n"
        )

print(f"VCF file created as '{output_filename}'")
