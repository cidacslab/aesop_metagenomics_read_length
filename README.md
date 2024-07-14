# Optimizing Next-Generation Sequencing Efficiency in Clinical Settings: Analysis of Read Length Impact on Cost and Performance

## Abstract

**Background:** The expansion of sequencing technologies as a result of the response to the COVID-19 pandemic enabled pathogen (meta)genomics to be deployed as a routine component of surveillance in many countries. Scaling genomic surveillance, however, come with associated costs in both equipment and sequencing reagents, which should be optimized. Here, we evaluate the cost efficiency and test performance of different read lengths in the identification of pathogens in metagenomic samples. We carefully evaluated performance metrics, financial implications, and time requirements relative to choices of 75 bp, 150 bp and 300 bp read lengths in pathogen identification metrics.  

**Results:** Our findings revealed that moving from 75 bp to 150 bp read length approximately doubles both the cost and sequencing time. Opting for 300 bp reads leads to four- and three-fold increases, respectively, in cost and sequencing time compared to 75 bp reads.  

**Conclusions:** During disease outbreaks situations, when responses are required, we suggest prioritizing shorter read lengths of 75 bp, as they provide satisfactory performance metrics while significantly reducing both costs (approximately two times lower) and, more importantly, -sequencing time (approximately three times faster). This practical approach enables projects to allocate their limited resources, streamline workflows, and maintain a response capability. 

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
