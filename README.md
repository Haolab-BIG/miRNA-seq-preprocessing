# miRNA-seq Pipeline

This unified miRNA-seq pipeline processes raw single-end FASTQ files through to differential expression results, using Singularity for reproducibility and supporting batch analysis of multiple samples.

## Workflow
<img width="1882" height="506" alt="CleanShot 2025-09-17 at 20 34 28@2x" src="https://github.com/user-attachments/assets/51d4db6a-87e7-4b8a-8224-84c8a5793eca" />


## Features

  * **Single Command Execution**: Executes the entire workflow—from FASTQ input and QC, through miRNA quantification, to differential expression analysis—with a single command.
  * **Reproducible**: All software (FastQC, Trim Galore, Bowtie1, miRDeep2, DESeq2) is encapsulated within a Singularity container (`miRNA.sif`), ensuring analysis is fully reproducible across different systems.
  * **Robust Quantification**: Leverages the well-established miRDeep2 package for accurate mapping and quantification of known miRNAs.
  * **Automated Reporting**: Generates a final, interactive MultiQC report summarizing quality control metrics across all samples and steps for easy assessment.

## Requirements

1.  **Recommended System Configuration**:

      * 8-core CPU
      * 32 GB RAM

2.  **Singularity**: Must be installed on your system. Below are detailed steps for installing on an Ubuntu 22.04 system. For other operating systems, please refer to the [official installation guide](https://www.google.com/search?q=https://docs.sylabs.io/guides/latest/user-guide/installation.html).

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
        	libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go (check for the latest version)
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Navigate to a suitable directory for downloading source code
        cd /tmp

        # Download the Singularity CE source code (check for the latest version)
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure, build, and install Singularity
        ./mconfig
        cd builddir
        make
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version
        ```

3.  **Pipeline Files**:

      * `run_pipeline.sh`
      * `miRNA.sif` (The Singularity container)

4.  **Reference Data**: A directory containing all necessary reference files.

## Setup

### 1\. Prepare the Sample Sheet

This is the most critical input file. Create a CSV file named `samplesheet.csv`. The pipeline is designed for **single-end** miRNA-seq data.

  * `sample`: A unique identifier for the sample (e.g., `Control_Rep1`). This name will be used for output subdirectories.
  * `condition`: The experimental group for the sample (e.g., `Control`, `Treated`). This is used for the DESeq2 design.
  * `fastq_path`: The **absolute path** to the single-end FASTQ file.

**Example `samplesheet.csv`:**

```csv
sample,condition,fastq_path
Control_Rep1,Control,/path/to/data/Control_Rep1.fastq.gz
Control_Rep2,Control,/path/to/data/Control_Rep2.fastq.gz
Treated_Rep1,Treated,/path/to/data/Treated_Rep1.fastq.gz
Treated_Rep2,Treated,/path/to/data/Treated_Rep2.fastq.gz
```

### 2\. Prepare the Reference Data

The pipeline requires several pre-built reference files for Human (hg38).

#### Create Reference Directory

Create a dedicated directory for all reference data:

```bash
mkdir -p reference_data
cd reference_data
```

#### Download Required Files

```bash
# 1. Download Genome FASTA (from GENCODE)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# 2. Download miRNA sequences (from miRBase)
# hairpin.fa contains precursor sequences
wget https://www.mirbase.org/download/hairpin.fa
# mature.fa contains mature miRNA sequences
wget https://www.mirbase.org/download/mature.fa

# Filter for human (hsa) miRNAs for efficiency
singularity exec ../miRNA.sif seqkit grep -r -p "^hsa" mature.fa -o mature_hsa.fa
singularity exec ../miRNA.sif seqkit grep -r -p "^hsa" hairpin.fa -o hairpin_hsa.fa
mv mature_hsa.fa mature.fa
mv hairpin_hsa.fa hairpin.fa

# 2. Download miRNA prediction files (from TargetScan)
wget https://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip
wget https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Info.default_predictions.txt.zip
gunzip miR_Family_Info.txt.zip
gunzip Predicted_Targets_Info.default_predictions.txt.zip
```


#### Build Bowtie1 Index

miRDeep2 uses Bowtie1 for mapping reads to the genome.

```bash
# The output prefix 'hsa' is used by the pipeline script
singularity exec ../miRNA.sif bowtie-build GRCh38.primary_assembly.genome.fa hsa
```

#### Final Reference Structure

Your `reference_data` directory should now look like this:

```
reference_data/
├── GRCh38.primary_assembly.genome.fa
├── hairpin.fa
├── mature.fa
├── miR_Family_Info.txt
├── Predicted_Targets_Info.default_predictions.txt
├── hsa.1.ebwt
├── hsa.2.ebwt
├── hsa.3.ebwt
├── hsa.4.ebwt
├── hsa.rev.1.ebwt
└── hsa.rev.2.ebwt
```

## Running

Execute the pipeline using a single command.

### Command Parameters

  * `-s`: Path to the sample sheet CSV file (required).
  * `-o`: Output directory path where results will be saved (required).
  * `-r`: Reference data directory (required).
  * `-c`: Path to the `miRNA.sif` Singularity container file (required).
  * `-L`: Control condition name for differential expression analysis (optional, if unset, DESeq2 uses alphabetical order).
  * `-t`: Number of threads to use for processing (optional, default: 8).

### Example Command

```bash
bash run_pipeline.sh \
  -s ./samplesheet.csv \
  -o ./mirna_project_results \
  -r ./reference_data \
  -c ./miRNA.sif \
  -L Control \
  -t 8
```

## Output Structure and Interpretation

After the pipeline completes successfully, the output directory (`mirna_project_results/`) will be organized as follows. The pipeline automatically cleans up most intermediate files, leaving only the key results.

```
./mirna_project_results/
├── Control_Rep1/
│   ├── Control_Rep1_mapped.fa
│   ├── Control_Rep1.arf
│   └── miRNAs_expressed.csv
├── Control_Rep2/
│   └── ... (same structure as Control_Rep1)
├── Treated_Rep1/
│   └── ...
├── Treated_Rep2/
│   └── ...
├── multiqc_report/
│   └── multiqc_report.html
├── deg_results.txt
├── miRNA_targets_genes.txt
└── merged_mirna_counts.csv
```

-----

### Aggregate Result Files

These files represent the final, combined analysis results from all samples.

   * **`deg_results.txt`**

      - **Content**: A tab-separated text file containing the results of the differential expression analysis from DESeq2. Each row corresponds to a gene, and columns typically include:
          - `baseMean`: Average normalized count across all samples.
          - `log2FoldChange`: The logarithm (base 2) of the fold change between the 'Treated' and 'Control' conditions. A positive value means the gene is upregulated in the 'Treated' group; a negative value means it is downregulated.
          - `lfcSE`: The standard error of the `log2FoldChange` estimate.
          - `stat`: The Wald statistic.
          - `pvalue`: The raw p-value for the statistical test.
          - `padj`: The p-value adjusted for multiple testing (e.g., using Benjamini-Hochberg correction).
      - **Application**: The deg_results.txt file is the primary input for many common downstream analyses and visualizations to interpret the biological meaning behind the differentially expressed genes like 'Volcano Plot'.
     <img width="1988" height="356" alt="CleanShot 2025-09-17 at 19 23 32@2x" src="https://github.com/user-attachments/assets/7b38d19d-2c83-4aab-8150-03e71dcf3942" />

  * **`merged_mirna_counts.csv`**

      * **Content**: A matrix of raw read counts. Rows are miRNAs, and columns are samples. This file is the direct input for the DESeq2 analysis.
      * **Application**: Can be used for custom quality control checks or visualizations, such as Principal Component Analysis (PCA) or sample-to-sample distance heatmaps to assess overall experiment consistency.

  * **`miRNA_targets_genes.txt`**

      * **Content**: A tab-separated text file containing the relationships between **miRNAs** and their **target genes**.  
           - **MiRBase_ID**: The standardized miRNA ID (e.g., *hsa-let-7g-5p*), usually from the **miRBase** database.  
           - **miR_Family**: The miRNA family information (e.g., *let-7-5p/98-5p*), used to reveal homologous or functionally related miRNA groups.  
           - **Gene_Symbol**: The official symbol of the predicted or validated target gene (e.g., *A1CF*, *ABCF1*), typically from a standard gene annotation database.  
      * **Application**: Provides a quick and simple list of the most important miRNAs from your analysis for reference or as input for other scripts.
      <img width="720" height="644" alt="CleanShot 2025-09-17 at 20 48 00@2x" src="https://github.com/user-attachments/assets/9126aac6-fefd-40be-bf6d-03bea1bebcc4" />

  * **`multiqc_report`**: Open `multiqc_report.html` in a web browser to explore all sections interactively.

  	- **Application**: This is the first file you should check to assess the overall quality of your sequencing data and the alignment process. It helps identify problematic samples (e.g., low alignment rate, high duplication) early on.

    	- **General Statistics**: A combined table summarizing important metrics for each sample:

    <img width="1924" height="770" alt="CleanShot 2025-09-17 at 19 19 05@2x" src="https://github.com/user-attachments/assets/4cf7ba35-669f-4247-b321-8f56e3ebba46" />

     - **FastQC**: Quality-control metrics on raw and trimmed reads, including  
      'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores',  
      'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content',  
      'Sequence Length Distribution', 'Sequence Duplication Levels',  
      'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content'.

     - **Sequence Quality Histograms**: The mean quality value across each base position in the read.  

       <img width="1914" height="1206" alt="CleanShot 2025-09-17 at 19 19 37@2x" src="https://github.com/user-attachments/assets/4f1955eb-ffb2-4a6b-929c-11400e7b239b" />


    	- **Botwie1**: Alignment statistics such as total reads, uniquely mapped reads, and multi-mapping rates:

     <img width="1910" height="950" alt="CleanShot 2025-09-17 at 19 21 11@2x" src="https://github.com/user-attachments/assets/877d517d-2e72-4037-b153-75e1ef8fbc46" />



-----

### Per-Sample Files (`Control_Rep1/`, etc.)

These directories contain key intermediate files for each individual sample, which are useful for debugging or more advanced analyses like novel miRNA prediction.

  * **`miRNAs_expressed.csv`**

      * **Content**: Per-sample quantification results from miRDeep2's `quantifier.pl` module, showing the raw read count for every known miRNA in that specific sample.
      * **Application**: Useful for inspecting the expression profile of a single sample. These individual files are parsed and combined to create the `merged_mirna_counts.csv` file.

  * **`*_mapped.fa`**

      * **Content**: A FASTA file containing all reads from the sample that successfully mapped to the reference genome.
      * **Application**: This file serves as a direct input for the quantification step. It can also be used, along with the `.arf` file, for novel miRNA prediction.

  * **`*.arf`**

      * **Content**: An Alignment Report Format (ARF) file that provides detailed coordinate information for each read's alignment in the `*_mapped.fa` file.
      * **Application**: Together, the `*_mapped.fa` and `*.arf` files are the required inputs for the main `miRDeep2.pl` script, which can be run manually to predict **novel miRNAs** that are not present in the current miRBase annotation.

## Video Tutorials


https://github.com/user-attachments/assets/8515662d-13f2-4c7d-9fda-d4bab8de2368

