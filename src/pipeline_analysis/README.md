# Execution of Analysis Pipeline

## Introduction

This section describes the steps followed to analyze the synthetic metagenomes, including removing ERCC reads, trimming adapters, filtering reads, annotating taxa, retrieving species-level taxa abundance, and normalizing the data.

## Steps

1. **Adapter Trimming and Quality Filtering**

    **Rationale:** Trimming adapters and filtering low-quality reads improve the accuracy of downstream analyses.

    **Command:**

    ```bash
    fastp -i input_ercc_removed.fastq -o output_trimmed.fastq
    ```

2. **Taxa Annotation**

    **Rationale:** Annotate reads to their respective taxa.

    **Command:**

    ```bash
    kraken2 --db path/to/kraken2_db --output kraken2_output.txt --report kraken2_report.txt --use-names output_human_removed.fastq
    ```

3. **Species-Level Taxa Abundance Retrieval**

    **Rationale:** Retrieve species-level taxa abundance for detailed analysis.

    **Command:**

    ```bash
    bracken -d path/to/kraken2_db -i kraken2_report.txt -o bracken_output.txt
    ```

4. **Calculate Taxa Annotation Confusion Matrix**

    **Rationale:** Calculate the confusion matrix of the taxa present in the sample compostition.

    **Command:**

    ```bash
    python calculate_confusion_matrix.py mock_sample_metadata.csv mock_sample.kreport mock_sample.kout mock_sample_confusion_matrix.csv
    ```