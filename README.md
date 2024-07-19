# Optimizing Next-Generation Sequencing Efficiency in Clinical Settings: Analysis of Read Length Impact on Cost and Performance

## Abstract

**Background**: The expansion of sequencing technologies as a result of the response to the COVID-19 pandemic enabled pathogen (meta)genomics to be deployed as a routine component of surveillance in many countries. Scaling genomic surveillance, however, comes with associated costs in both equipment and sequencing reagents, which should be optimized. Here, we evaluate the cost efficiency and performance of different read lengths in identifying pathogens in metagenomic samples. We carefully evaluated performance metrics, costs, and time requirements relative to choices of 75 bp, 150 bp and 300 bp read lengths in pathogen identification. 

**Results**: Our findings revealed that moving from 75 bp to 150 bp read length approximately doubles both the cost and sequencing time. Opting for 300 bp reads leads to four- and three-fold increases, respectively, in cost and sequencing time compared to 75 bp reads. For viral pathogen detection, the sensitivity median ranged from 97.9% with 75 bp reads to 100% with 150 or 300 bp reads. However, bacterial pathogens detection was less effective with shorter reads: 76% with 75 bp, 90% with 150 bp, and 94.3% with 300 bp reads. These findings were consistent across different levels of taxa abundance.

**Conclusions**: During disease outbreak situations, when swift responses are required for pathogen identification, we suggest prioritizing 75 bp read lengths. Shorter reads enable quicker sequencing times (approximately three times faster) and reduce costs (approximately two times lower). Despite the shorter read length, the performance in terms of precision is comparable to that of longer reads across most viral and bacterial taxa, while sensitivity can be more variable, especially if bacterial identification is aimed. This practical approach allows better use of resources, enabling the sequencing of more samples using streamlined workflows, while maintaining a reliable response capability.


## Methods

Our work performed the following steps:

1. [**Generation of Synthetic Metagenomes**](/src/metagenome_creation/README.md)
    1. Defining the metagenome composition
    2. Defining each taxon abundance
    3. Collecting the synthetic sample taxonomic data
    4. Downloading the genomes of these taxa
    5. Generation of the synthetic metagenomes
2. [**Execution of Analysis Pipeline**](/src/pipeline_analysis/README.md)
    1. Adapter Trimming and Quality Filtering
    2. Taxa Annotation
    3. Species-Level Taxa Abundance Retrieval
    4. Calculate Each Taxa Confusion Matrix
3. [**Creating and Plotting Results**](/src/paper_figures/README.md)


## Installation

Install the necessary software using the following commands:

```bash
# Install Fastp
conda install -c bioconda fastp

# Install Kraken2
conda install -c bioconda kraken2

# Install Bracken
conda install -c bioconda bracken
```

## Usage

1. **Clone the repository**

```bash
git clone https://github.com/your_username/your_repository.git
cd your_repository
```

2. **Execute the pipeline**

Follow the steps detailed in our [METHODS](#methods)


## Citation

If you use this pipeline in your research, please cite the following paper:


> Meirelles, P. M.; Viana, P. A. B.; Tschoeke, D. A.; de Moraes, L.; Amorim, L.; Barral-Netto, M.; Khouri, R.; Ramos, P. I. P. (2024). Optimizing Next-Generation Sequencing Efficiency in Clinical Settings: Analysis of Read Length Impact on Cost and Performance.

* Corresponding Author: Pedro M Meirelles (pmeirelles@ufba.br)
* On any code issues, correspond to: Pablo Viana (pablo.alessandro@gmail.com)
