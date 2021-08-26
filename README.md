# Pomfrey
This branch is for running with [Pgx_module](https://github.com/JoelAAs/pgx_module/tree/pgx_pomfrey) on branch pgx_pomfrey.
##### Hematology Twist Pipeline

To run the pipeline you need Snakemake and Singularitys installed. At Uppsala it is used together with slurm-drmaa to submit on the local HPC. If Horizon Myeloid DNA Reference Standard is used it should be named HD829 to be processed separately and not hold up the pipeline.


### Installation
```
git clone https://github.com/clinical-genomics-uppsala/pomfrey
cd pomfrey 
git checkout pomfey_pgx
git clone https://github.com/joelaas/pgx_module
cd pgx_module 
git checkout pgx_pomfrey

cd envs
chmod +x get_cointainers.sh
sudo ./get_containers.sh
```

For pgx-specific parameters please read: [README](https://github.com/JoelAAs/pgx_module/tree/pgx_pomfrey#readme)

### Files & Caches
- **SampleSheet.csv**: A csv-file with produced when demultiplexing on Illumina machine. Is used to order samples in MultiQC table. The script only use the column *Sample_Name*. Lines needed are in the file except the actual sample-lines are:
    ```sh
    [Data]
    Sample_ID,Sample_Name,Description,index,I7_Index_ID,index2,I5_Index_ID,Sample_Project
    ```
- **Bedfile**: Four columns: chr, start, stop, regionname. No header. The bed-file is used as both target and bait intervals in picard HsMetric for bases on target stats. The fourth column is used in CARTools to identify low coverage regions.
- **Pindel bedfile**: Since pindel is quite slow a smaller bedfile with limited regions is needed to run pindel.
- **Intervals-file**: Corresponding to the bedfile for picard jobs. Can be generated in GATK4 with:
    `gatk fileBedToIntervalList --INPUT $bedfile --OUTPUT $interval_output --SEQUENCE_DICTIONARY $reference_dict`
- **Reference**: Fasta refernce index both with bwa index and .fai.
- **Artefact file**: Known artefacts in bedformat
- **Germline file**:
- **Pindel artefact file**: Known pindel artefacts, chr, pos, gene.
- **COSMIC hemato counts**: COSMIC file with number of hemato hits. Columns (tab seperated): GENE_NAME, ACCESSION_NUMBER, GENE_CDS_LENGTH, HGNC_ID PRIMARY_SITE, PRIMARY_HISTOLOGY,  GENOMIC_MUTATION_ID, LEGACY_MUTATION_ID, MUTATION_ID, MUTATION_CDS, MUTATION_AA, MUTATION_GENOME_POSITION, MUTATION_STRAND, SNP, MUTATION_SOMATIC_STATUS, N_observations, Chromosome, chr, Start, End, Build
- **Hotspot list**: List of regions where higher coverage is important. Each region and its coverage is listed in the sheet "HotSpot" in the xlsx-file
- **VEP cache**: Need to download cache for vep to run. Read more on the different versions of [vep-cache](https://m.ensembl.org/info/docs/tools/vep/script/vep_cache.html).
    `singularity exec --bind $PWD vep-container.simg perl /opt/vep/src/ensembl-vep/INSTALL.pl -s homo_sapiens_refseq --CACHEDIR vep-data-99.0/ -a c --ASSEMBLY GRCh37`


### Singularity Containers
| Program | Version | Source |
| ------- | ------- | ------ |
| Bcftools |1.9 | docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_8	|
| Cutadapt | 2.5 |	docker://quay.io/biocontainers/cutadapt:2.5--py37h516909a_0 |
| Fastqc | 0.11.8 | docker://quay.io/biocontainers/fastqc:0.11.8--1	|
| Freebayes | 1.3.1 |docker://quay.io/biocontainers/freebayes:1.3.1--py37h56106d0_0 |
| GATK4 | 4.1.7.0 | docker://broadinstitute/gatk:4.1.7.0 |
| MultiQC| 1.7 |	docker://quay.io/biocontainers/multiqc:1.7--py_3	|
| Pindel | 0.2.58	| docker://shuangbroad/pindel:v0.2.5b8 |
| Vardict-java | 1.7.0 | docker://quay.io/biocontainers/vardict-java:1.7.0--0	|
| Vep | 99 |docker://ensemblorg/ensembl-vep:release_99.0	|
| Vt | 0.57721 | docker://quay.io/biocontainers/vt:0.57721--hdf88d34_2	|

|Program| Comment| Source|
| ----- | ------ | ----- |
| bcbio-ensemble-recall | |		bcbio-variation-recall.simg |
| Bwa & Samtools & Picard|  v.7.17 &  v.1.9 & v.2.20.1 | bwa-snakemakewrapper.def |
| CARTool |In-house program avaialble at [Github](https://github.com/anod6351/CARtool) | CARTools-200206.simg |
| Igv v.2.4.10 & xvfb | Not working to load hg19 genome into container | igv.def |
| Pisces v.5.2.11.163 |Microsoft dotnet v.2.1| Pisces_singularity.def	|
| Python3 | Libraries: csv, pysam, xlsxwriter,date, itemgetter, yaml|	python-pysam.def |

### Config Files
#### Samples Config
Yaml config file with samples to process and local variables and files for the pipeline. Configfile name has to be ${sequencerun}_config.yaml
```
programdir:
  dir: "${PATH_TO_POMFREY}"

reference:
    ref: "" #Path to fasta ref and .fai
    bwa: "" #Path to bwa indexed ref

configCache:
    multiqc: "${PATH_TO_POMFREY}/src/report/multiqc_config.yaml"
    vep: "" #Path to downloaded VEP cache
    hemato: "" #Path to COSMIC file for number of hematology hits. Downloaded from COSMIC site..
    variantlist: "" #Path to file where all variants are written to later create artefact and germlinefilter files.

bed:
    bedfile: "" #Path to main bedfile
    intervals: "" Path to interval file corresponding to main bedfile
    pindel: ""  #Path to pindel bedfile
    indelartefact: "" # Path to indel artefacts since all indels from vardict and mutect2 are included
    pindelArtefact: ""  #Path to pindel artefact file
    cartool: "" #Path to (main) bedfile or different if interested in the coverage of different regions.
    hotspot: "" #Path to hotpspotlist
    artefact: "" #Path to artefact filter file
    germline: "" #Path to germline filter file

CNV:
    PoN: "" # Path to you readCountPoN.hdf5 file for the panel of normals
    bedPoN: "" # Bedfile used for PoN, could be normal bedfile
    interval: "" #interval list
    cyto: "" # Translation from cyto coordinate to genomic position [chr start end cytoCoord]

singularitys:
    cutadapt: ""
    bwa: ""
    fastqc: ""
    cartool: ""
    bcftools: ""
    freebayes: ""
    pisces: ""
    vardict: ""
    gatk4: ""
    pindel: ""
    vep: ""
    recall: ""
    python: ""
    vt: ""
    igv: ""
    multiqc: ""

cartool:
    cov: "100 200 1000" #Coverage limits, first number minCov, second for hotspotlist, third wishful

methods:   # The (trust) order of vcfs into ensemble recall
    mutect2: "mutect2"
    vardict: "vardict"
    pisces: "pisces"
    freebayes: "freebayes"

seqID:
    sequencerun: ""  

samples:
   "${sampleName}": "${path_to_R1_fastq_gz}"

```

#### Cluster Config
Json file with config for submission on HPC. Need to be specified to suit you HPC. See cluster-config.json for example.
### Snakemake command
`
snakemake -p -j ${max_nr_jobs_submitted} --drmaa "-A ${project} -p core -t {cluster.time} -n {cluster.n} --nodes=1-1 " --use-singularity --cluster-config ${cluster_config} -s ${PATH_TO_POMFREY}/src/somaticPipeline.smk --singularity-args " --cleanenv --bind /data/ --bind /projects/ " --configfile ${sample_config}
`
