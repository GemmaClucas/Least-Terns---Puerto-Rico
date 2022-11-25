MiFish bioinformatics
================
Gemma Clucas
2022-11-25

## R Markdown

## 1. Import the data into Qiime2

Load qiime environment and cd to correct directory.

    cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/Least-Terns---Puerto-Rico
    conda activate qiime2-2021.4

The data is saved on my solid state hard drive, but I forgot to bring it
home, so I redownloaded it to this directory. I should delete these as
soon as I’ve got the data into a qiime archive.

    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Users/gemmaclucas/GitHub/Fecal_metabarcoding/Least-Terns---Puerto-Rico/reads/ \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path demux_plate1.qza

## 2. Trim primers with cutadapt plugin

F primer: GTCGGTAAAACTCGTGCCAGC (21 bp)  
R primer: CATAGTGGGGTATCTAATCCCAGTTTG (27 bp)

### Trim 3’ ends first

At the 3’ end of the read, the primer will have been read through after
reading the MiFish sequence. I need to be looking for the reverse
complement of the reverse primer in R1 (`—p-adapter-f`) and the reverse
complement of the forward primer in R2 (`—p-adapter-r`)

F primer reverse complement: GCTGGCACGAGTTTTACCGAC  
R primer reverse complement: CAAACTGGGATTAGATACCCCACTATG

    for K in {1..1}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux_Plate$K.qza \
        --p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
        --p-adapter-r GCTGGCACGAGTTTTACCGAC \
        --o-trimmed-sequences trimd_Plate$K.qza \
        --verbose > cutadapt_out_Plate$K.txt
    done

To see how much data passed the filter for each sample:

    grep "Total written (filtered):" cutadapt_out_Plate1.txt 

Between 50 - 77% is passing the filters - should be good

### Trim 5’ ends of reads

All R1 should begin with the forward primer: GTCGGTAAAACTCGTGCCAGC (21
bases).  
All R2 should begin with the reverse primer: CATAGTGGGGTATCTAATCCCAGTTTG
(27 bases).

Trim these with the following commands:

    for K in {1..1}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences trimd_Plate$K.qza \
        --p-front-f GTCGGTAAAACTCGTGCCAGC \
        --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
        --o-trimmed-sequences trimd2_Plate$K.qza \
        --verbose > cutadapt_out2_Plate$K.txt
    done

About 87% of the data is passing the filter very consistently this time,
which is the usual amount.

    grep "Total written (filtered):" cutadapt_out2_Plate1.txt

## 3. Denoise with dada2

I am going to use the same settings that I used for the 2017 and 2018
tern fecal samples here, except I need to add the `--p-min-overlap`
parameter, otherwise I seem to be getting a load of rubbish reads which
are 250bp long and start with long strings of Cs. This is only available
in qiime2-2021.4 and later. I initially tried specifying an overlap of
30, but that didn’t seem like enough as I was still getting junk
sequences, but I think 50 is working well now.

Note, this step is pretty slow to run, a whole plate takes about 40 mins
(but that changes depending on sequencing depth).

    for K in {1..1}; do
      qiime dada2 denoise-paired \
        --i-demultiplexed-seqs trimd2_Plate$K.qza \
        --p-trunc-len-f 133 \
        --p-trunc-len-r 138 \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-min-overlap 50 \
        --p-n-threads 12 \
        --o-representative-sequences rep-seqs_Plate$K \
        --o-table table_Plate$K \
        --o-denoising-stats denoise_Plate$K
    done

Create visualizations for the denoising stats.

    for K in {1..1}; do  
      qiime metadata tabulate\
        --m-input-file denoise_Plate$K.qza\
        --o-visualization denoise_Plate$K.qzv
    done

Looks like sequencing went really well - most samples had more than 90%
of sequences retained.

## 4. Assign taxonomy

I’m going to use the database I created for Will’s ATPU project, which
contains all fish 12S sequences on GenBank in 2021, and Devin’s blast
method. The blast method uses an older version of qiime.

    conda activate qiime2-2019.4

    ./mktaxa.py ncbi-refseqs-withHuman.qza \
      ncbi-taxonomy-withHuman.qza \
      rep-seqs_Plate1.qza
      
    qiime metadata tabulate \
      --m-input-file superblast_taxonomy.qza \
      --o-visualization superblast_taxonomy

## 5. Make barplots

MAKE METADATA FILE

    qiime taxa barplot \
      --i-table table_merged.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_before_filtering.qzv