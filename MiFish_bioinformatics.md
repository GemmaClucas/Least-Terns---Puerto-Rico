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
      
    qiime feature-table tabulate-seqs \
      --i-data rep-seqs_Plate1.qza \
      --o-visualization rep-seqs_Plate1

## 5. Make barplots

    conda activate qiime2-2021.4

    qiime taxa barplot \
      --i-table table_Plate1.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_before_filtering.qzv

Not too much bird or human DNA in here, but lots of fish diversity. It
will take me a while to sort through all these IDs.

E.g. Thysochromis ansorgii is the most common species but it is only a
90% match. Seems like this is possibly a fish in the Antherinidae
family, the old word silversides, which are found in coastal waters in
tropical and temperate latitides. However, these seem to be mostly found
in the Old World, so probably safest to just say it is something in the
order Antheriniformes (silversides).

## 6. Remove non-food reads

I need to filter out any sequences from the bird, mammals, and
unnassigned sequences before rarefying.

    qiime taxa filter-table \
      --i-table table_Plate1.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --p-exclude Unassigned,Aves,Mammalia \
      --o-filtered-table table_Plate1_noBirdsMammalsUnassigned.qza
      
    qiime feature-table summarize \
        --i-table table_Plate1_noBirdsMammalsUnassigned.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization table_Plate1_noBirdsMammalsUnassigned

## 7. Calculate rarefaction curves

Samples range from 3 reads to 500,000. I don’t want to go lower than 100
and I’ll set the upper limit to 20,000.

    qiime taxa collapse \
      --i-table table_Plate1_noBirdsMammalsUnassigned.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --p-level 7 \
      --o-collapsed-table table_Plate1_noBirdsMammalsUnassigned_collapsed.qza

    qiime diversity alpha-rarefaction \
      --i-table table_Plate1_noBirdsMammalsUnassigned_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 100 \
      --p-max-depth 20000 \
      --o-visualization alpha-rarefaction-100-20000

It does look like there is high diversity in these samples so I probably
need a higher depth than in the Gulf of Maine terns. If you look at all
samples together, then they plateau at about 6500 sequences -\> redo the
rarefaction from 100 to 10,000

    qiime diversity alpha-rarefaction \
      --i-table table_Plate1_noBirdsMammalsUnassigned_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 100 \
      --p-max-depth 8000 \
      --o-visualization alpha-rarefaction-100-8000

By 5500 samples have plateaued.

## 8. Rarefy to a depth of 5500

And redo the barplots.

    qiime feature-table rarefy \
      --i-table table_Plate1_noBirdsMammalsUnassigned.qza \
      --p-sampling-depth 5500 \
      --o-rarefied-table table_Plate1_noBirdsMammalsUnassigned_rarefied5500 
      
    qiime taxa barplot\
      --i-table table_Plate1_noBirdsMammalsUnassigned_rarefied5500.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_Plate1_noBirdsMammalsUnassigned_rarefied5500

## 9. Check taxonomy assignments

Had to make a lot of edits as it seems like there are a lot of sequences
that are good matches to many species.

-   Harengula jaguana - good, sacled herring.
-   Thysochromis ansorgii - changed to order Atheriniformes sp.,
    silversides.
-   Jenkinsia lamprotaenia - good, dwarf round herring.
-   Ctenogobius boleosoma - changed to subfamily Gobionellinae, gobies.
-   Kyphosus incisor - changed to Kyphosus sp., sea chubs.
-   Hyporhamphus yuri - changed to Hemiramphidae, halfbeaks (90% is
    closest match).
-   Elagatis bipinnulata - good, rainbow runner.
-   Cheilopogon antoncichi - changed to family Exocoetidae, flying fish.
-   Carangoides bartholomaei - changed to Caranx crysos, blue runner.
-   Anchoa;s\_\_hepsetus - changed to Anchoa sp., anchovy.
-   Seriola rivoliana - good, longfin yellowtail.
-   Cephalopholis;s\_\_leopardus - changed to Eucinostomus gula, Jenny
    mojarra
-   Hemiramphus lutkei - changed to Hemiramphus sp., halfbeaks (99%
    match to multiple species in the genus)
-   Kyphosus bigibbus - added Kyphosus sp., sea chubs.
-   Mulloidichthys vanicolensis - changed to Mulloidichthys sp.,
    goatfish.
-   Stegastes planifrons - changed to Pomacentridae sp., damselfishes.
-   Coryphaena hippurus - good, common dolphonfish.
-   Craterocephalus;s\_\_fluviatilis - added to the Atherininae,
    silversides.
-   Engraulis encrasicolus - chnaged to Engraulis sp., anchovy.
-   Eleutheronema rhadinum - changed to Polynemidae sp., threadfins.
-   Tylosurus crocodilus crocodilus - good, hound needlefish.
-   Hemiramphus far - added to Hemiramphus sp. as not found in Atlantic.
-   Abudefduf saxatilis - changed to Abudefduf sp., seargent majors
-   Anchoa delicatissima - added to Anchoa sp., anchovy
-   Exocoetus monocirrhus - added to family Exocoetidae, flying fish.
-   Trachinotus;s\_\_mookalee - changed to Trachinotus sp., pompanos.
-   Hyporhamphus quoyi - added to Hemiramphidae, halfbeaks.
-   Hirundichthys oxycephalus - added to family Exocoetidae, flying
    fish.
-   Gerres cinereus - good, yellow-fin mojarra.
-   Cheilopogon doederleinii - added to family Exocoetidae, flying fish.
-   Gambusia affinis - changed to Gambusia sp., mosquitofishes.
-   Exocoetus volitans - added to family Exocoetidae, flying fish.

I sent the data over to Luis as both presence/absence and RRA
