import pandas as pd
from sklearn.metrics import confusion_matrix

def load_kraken_output(kraken_file):
    # Load Kraken output file
    kraken_data = pd.read_csv(kraken_file, sep='\t', header=None, names=['classified', 'seq_id', 'tax_id'])
    kraken_data['seq_id'] = kraken_data['seq_id'].str.split('|').str[1]  # Extract sequence ID part
    return kraken_data

def load_metadata(metadata_file):
    # Load metadata file mapping sequence IDs to actual tax IDs
    metadata = pd.read_csv(metadata_file, sep='\t')
    return metadata

def merge_data(kraken_data, metadata):
    # Merge Kraken data with metadata to get the actual tax IDs
    merged_data = pd.merge(kraken_data, metadata, left_on='seq_id', right_on='accession', how='inner')
    return merged_data[['seq_id', 'tax_id_x', 'tax_id_y']].rename(columns={'tax_id_x': 'kraken_tax_id', 'tax_id_y': 'true_tax_id'})

def calculate_confusion_matrix(merged_data):
    # Calculate confusion matrix
    y_true = merged_data['true_tax_id']
    y_pred = merged_data['kraken_tax_id']
    labels = sorted(set(y_true) | set(y_pred))
    conf_matrix = confusion_matrix(y_true, y_pred, labels=labels)
    return conf_matrix, labels

def save_confusion_matrix(conf_matrix, labels, output_file):
    # Save confusion matrix to a CSV file
    conf_matrix_df = pd.DataFrame(conf_matrix, index=labels, columns=labels)
    conf_matrix_df.to_csv(output_file)

def main():
    kraken_file = 'kraken_output.txt'
    metadata_file = 'metadata.csv'
    output_file = 'confusion_matrix.csv'
    
    # Load data
    kraken_data = load_kraken_output(kraken_file)
    metadata = load_metadata(metadata_file)
    
    # Merge data and calculate confusion matrix
    merged_data = merge_data(kraken_data, metadata)
    conf_matrix, labels = calculate_confusion_matrix(merged_data)
    
    # Save confusion matrix
    save_confusion_matrix(conf_matrix, labels, output_file)
    print(f'Confusion matrix saved to {output_file}')

if __name__ == "__main__":
    main()
