# Digester-MultiSequencing

**Description:**
Repository for scripts and resources used for genome-centric metagenomics of anaerobic digester sludge (and Zymo Mock) with different sequencing approaches (Nanopore R9/R10, PacBio HiFi, Illumina short reads).
<br/>

**Quick links:**
* Pre-print for the paper can be accessed [here](https://www.biorxiv.org/content/10.1101/2021.10.27.466057v2).
* Sequencing read datasets are available at ENA for the [Zymo](https://www.ebi.ac.uk/ena/browser/view/PRJEB48692) and [anaerobic digester](https://www.ebi.ac.uk/ena/browser/view/PRJEB48021) samples.
* Anaerobic digester MAGs and Zymo assembly sequences can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.17008801). 
<br/>

# Summary of using Nanopore R10.4 for genome-centric metagenomics

* A high complexity metagenomic sample (anaerobic digester sludge) was sequenced with Nanopore R10.4 as well as Illumina Miseq, Nanopore R9.4.1 and PacBio HiFi to compare the different sequencing platforms. Overview of bionformatic processing steps is presented below:

<img src="https://github.com/Serka-M/Digester-MultiSequencing/blob/main/code/figs/mags_r104_workflow.svg" alt="AD-bioinf" style="zoom:100%;" />
<br/>

* Using PacBio HiFi assembly polished with Illumina reads as a reference, Nanopore R10.4 assembly was found to feature improved homopolymer calling, compared to Nanopore R9.4.1, especially for guanines and cytosines:

<img src="https://github.com/Serka-M/Digester-MultiSequencing/blob/main/code/figs/mags_r104_hp.svg" alt="AD-bioinf" style="zoom:100%;" />
<br/>

* Improvement of homopolymer calling for Nanopore R10.4 data is significant for metagenomics as most microbial genomes do not feature many homopolymers above the length of 10. To illustrate this, homopolymer rates were counted in genomes from RefSeq database:

<img src="https://github.com/Serka-M/Digester-MultiSequencing/blob/main/code/figs/hp_refseq.svg" alt="AD-bioinf" style="zoom:100%;" />
<br/>
