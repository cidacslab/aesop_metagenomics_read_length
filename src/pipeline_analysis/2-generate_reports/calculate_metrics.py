import pandas as pd

# Load the necessary data
ground_truth_hierarchy_path = 'results/mocks_throat_based/metadata_tree/throat_non_pathogen_02.csv'
classification_data_path = 'results/mocks_throat_based/performance_metrics/throat_non_pathogen_02_300bp_reads_out.csv'
kreport_path = 'results/mocks_throat_based/pipeline_outputs/2-kraken_results/throat_non_pathogen_02_300bp_reads.kreport'
final_output_path = 'results/mocks_throat_based/performance_metrics/throat_non_pathogen_02_300bp_reads_metrics.csv'

ground_truth_data = pd.read_csv(ground_truth_hierarchy_path)

classification_data = pd.read_csv(classification_data_path)
classification_data.columns = ['accession', 'kraken_taxid', 'read_count']

kreport_data = pd.read_csv(kreport_path, delimiter='\t', header=None, names=[
    'percent_reads', 'num_reads', 'num_reads_this_level', 'taxonomic_rank', 'NCBI_taxonomic_ID', 'scientific_name'
])

# Merge ground truth with classification data to include total reads
total_reads = classification_data.groupby('accession')['read_count'].sum().reset_index()
total_reads.columns = ['accession', 'total_accession_reads']

ground_truth_data = ground_truth_data.merge(total_reads, on='accession', how='left')

# Extract genus and species taxonomic IDs
ground_truth_data['genus_name'] = ground_truth_data['genus'].apply(lambda x: eval(x)[0] if pd.notna(x) else None).astype(str)
ground_truth_data['genus_taxid'] = ground_truth_data['genus'].apply(lambda x: int(eval(x)[1]) if pd.notna(x) else None).astype(str)
ground_truth_data['species_name'] = ground_truth_data['species'].apply(lambda x: eval(x)[0] if pd.notna(x) else None).astype(str)
ground_truth_data['species_taxid'] = ground_truth_data['species'].apply(lambda x: int(eval(x)[1]) if pd.notna(x) else None).astype(str)

# Create confusion matrices from the kreport data
classification_data['kraken_taxid'] = classification_data['kraken_taxid'].astype(int)

# Function to calculate TP, FP, FN, and TN for each accession
def calculate_accession_metrics(accession_id, taxid, rank):
    # Filter classification data for the current accession
    accession_data = classification_data[classification_data['accession'] == accession_id]
    
    # Create a confusion matrix for the current accession
    confusion_matrix = pd.crosstab(
        accession_data['kraken_taxid'], 
        [accession_data[rank]],
        values=accession_data['read_count'],
        aggfunc='sum', 
        dropna=False
    ).fillna(0).astype(int)
    
    TP = int(confusion_matrix.loc[taxid, taxid]) if taxid in confusion_matrix.index and taxid in confusion_matrix.columns else 0
    FP = int(confusion_matrix[taxid].sum() - TP) if taxid in confusion_matrix.columns else 0
    FN = int(confusion_matrix.loc[taxid].sum() - TP) if taxid in confusion_matrix.index else 0
    TN = int(confusion_matrix.values.sum() - (TP + FP + FN))
    
    return {'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN}

# Calculate metrics for each accession
results = []

for idx, row in ground_truth_data.iterrows():
    accession_id = row['accession']
    
    genus_taxid = row['genus_taxid']
    species_taxid = row['species_taxid']
    
    genus_metrics = calculate_accession_metrics(accession_id, genus_taxid, 'genus_taxid')
    species_metrics = calculate_accession_metrics(accession_id, species_taxid, 'species_taxid')
    
    results.append({
        'accession_id': accession_id,
        'total_accession_reads': row['total_accession_reads'],
        'genus_name': row['genus_name'],
        'genus_taxid': genus_taxid,
        'genus_true_positive': genus_metrics['TP'],
        'genus_true_negative': genus_metrics['TN'],
        'genus_false_positive': genus_metrics['FP'],
        'genus_false_negative': genus_metrics['FN'],
        'species_name': row['species_name'],
        'species_taxid': species_taxid,
        'species_true_positive': species_metrics['TP'],
        'species_true_negative': species_metrics['TN'],
        'species_false_positive': species_metrics['FP'],
        'species_false_negative': species_metrics['FN']
    })

# Convert results to DataFrame
final_data = pd.DataFrame(results)

# Select and reorder the required columns
final_columns = [
    'accession_id', 'total_accession_reads', 'genus_name', 'genus_taxid', 'genus_true_positive', 'genus_true_negative',
    'genus_false_positive', 'genus_false_negative', 'species_name', 'species_taxid', 'species_true_positive',
    'species_true_negative', 'species_false_positive', 'species_false_negative'
]
final_data = final_data[final_columns]

# Save the final output to CSV
final_data.to_csv(final_output_path, index=False)

print("Final metrics output saved to:", final_output_path)
print(final_data.head())
