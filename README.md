# Lynch-Sybdrome-Detection
A bash script using fastq file to detect Lynch syndrome (hereditary nonpolyposis colorectal cancer (HNPCC)).
This project aims to detect Lynch syndrome.

## Data Files

The fastq files used for analysis can be found in the following location:

- Description: [Link to fastq files NCBI Database](https://www.ncbi.nlm.nih.gov/sra)

## Reference Genome

The reference genome files required for indexing can be found in the following location:

- Path: `/path/to/reference2/folder`
- Description: [Link to reference genome files](https://drive.google.com/file/d/1-85dOBkNpqRyP99yt8P-IO4p7wFnRPJ_/view?usp=sharing)

## Tools

The required tools for the analysis are available in the following location:

- Path: `/path/to/tools/folder`
- Description: [Link to tools](https://drive.google.com/drive/folders/18q7gafdHOHe7Hoo-ZVWeXwOExroNmBbN?usp=sharing)

## Analysis Pipeline

The analysis pipeline involves the following steps:

1. Downloading fastq files from the provided link.
2. Checking the quality of the fastq files using `fastqc`.
3. Downloading and setting up Trimmomatic for trimming the reads.
4. Indexing the human genome using `bwa index`.
5. Aligning the reads to the reference genome using `bwa mem`.
6. Converting the alignment to BAM format using `samtools view`.
7. Sorting and indexing the BAM file using `samtools sort` and `samtools index`.
8. Calling variants using `bcftools mpileup` and `bcftools call`.
9. Performing annotation using `convert2annovar.pl` and `annotate_variation.pl`.
10. Filtering specific positions using `grep`.

## Usage

To run the analysis, follow these steps:

1. Clone the repository:

   ```bash
   git clone https://github.com/mariamadel41/Lynch-Sybdrome-Detection.git

2. Set up the necessary data files, reference genome, and tools as described above.

3. Modify the script files if required to match your specific file paths and configurations.

4. Run the main script:

   ```bash

    bash run_analysis.sh /path/file.fastq

This will execute the entire analysis pipeline.

## License

This project is licensed under the MIT License.

Feel free to customize the content and formatting of the `README.md` file according to your project's specific requirements.

Please make sure to replace the placeholder paths and descriptions with the actual ones relevant to your project.
