# Harper_et_al_2018

Data processing workflow and supplementary data for Harper *et al.* 2018:

Needle in a haystack? A comparison of eDNA metabarcoding and targeted qPCR for detection of great crested newt (*Triturus cristatus*)

and

Abiotic and biotic determinants of great crested newt at the pondscape using environmental DNA


## Contents

Notebooks to create curated reference databases used in analyses (databases also available in Genbank format) [(here)](https://github.com/HullUni-bioinformatics/Harper_et_al_2018/tree/master/Reference%20databases)

SRA accession numbers for raw Illumina data (here)

Taxonomic assignment results [(here)](https://github.com/HullUni-bioinformatics/Harper_et_al_2018/tree/master/Data)

R scripts used to analyse metaBEAT output and produce figures (here)



## Instructions to set up dependencies for data processing and analyses

To facilitate full reproducibility of our analyses, we provide Jupyter notebooks illustrating our workflow and all necessary supplementary data in this repository.

Illumina data was processed (from raw reads to taxonomic assignment) using the metaBEAT pipeline. The pipeline relies on a range of open bioinformatics tools, which we have wrapped up in a self contained docker image which includes all necessary dependencies here.



## Setting up the environment

In order to retrieve supplementary data (reference sequences etc.), start by cloning this repository to your current directory:

git clone --recursive https://github.com/HullUni-bioinformatics/Harper_et_al_2018.git

In order to make use of our self contained analysis environment, you will have to install Docker on your computer. Docker is compatible with all major operating systems. See the Docker documenation for details. On Ubuntu installing Docker should be as easy as:

sudo apt-get install docker.io

Once Docker is installed you can enter the environment by typing:

sudo docker run -i -t --net=host --name metaBEAT -v $(pwd):/home/working chrishah/metabeat /bin/bash

This will download the metaBEAT image (if it's not yet present on your computer) and enter the 'container' i.e. the self contained environment (NB: sudo may be necessary in some cases). With the above command the container's directory /home/working will be mounted to your current working directory (as instructed by $(pwd)). In other words, anything you do in the container's /home/working directory will be synced with your current working directory on your local machine.



## Data processing workflow as Jupyter notebooks

Raw illumina data has been deposited on the NCBI SRA (BioProject: ...; BioSample accession: ...; Sequence Read Archive accessions: ...). The sample specific accessions can be found here. Before following the workflow for data processing, you'll need to download the raw reads from SRA. To download the raw read data you can follow the steps in this Jupyter notebook.

With the data in place, you should be able to fully reproduce our analyses by following the steps outlined in the Jupyter notebooks.

The workflow illustrated in the notebooks assumes that the raw Illumina data is present in a directory raw_reads at the base of the repository structure and that the files are named according to the following convention: 'sampleID-marker', followed by '_R1' or '_R2' to identify the forward/reverse read file respectively. sampleID must corresponds to the first column in the file Sample_accessions.tsv here.
