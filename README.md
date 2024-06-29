# New virus detection based on the optimal natural metric

## Abstract
During the COVID-19 pandemic, the highly variable SARS-CoV-2 virus frequently undergoes mutations, leading to the emergence of new variants that present novel threats to public health. The determination of these variants often relies on manual definition based on local sequence characteristics, resulting in delays in their detection relative to their actual emergence. In this study, we propose an algorithm for the automatic identification of novel variants. By leveraging the optimal natural metric for viruses based on an alignment-free perspective to measure distances between sequences, we devise a hypothesis testing framework to determine whether a given viral sequence belongs to a novel variant. Our method demonstrates high accuracy, achieving nearly 100% precision not only in identifying new variants of SARS-CoV-2 and HIV-1, but also in detecting novel genus in Orthocoronavirinae. This approach holds promise for timely surveillance and management of emerging viral threats in the field of public health.

## Environment Setup
- Python 3
- Numpy
- Biopython

## Data in this study
The csv file containing sequence information and the fasta file (zip format )containing sequences are all available.

## Generate distance matrix under the optimal metric

Let Fasta_Filename = the name of your fasta file, then run the code "metric.py".

## Detect new variants

Let Csv_Filename = the name of your file containing variant labels (npy format), then run the code "algorithm1.py".
