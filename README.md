# Blowfly project

The goal of the blowfly projetc is 1. update a newly chromosomal scale genome assembly 2. Using the genome as reference 

# Blowfly Genome Assembly

The script was written by Sheng-Hao (Mayson). The pipeline shown is to generate a contiguous and comprehensive chromosomal scale level of Hi-C intergated genome assembly, providing valuable insights into the genomic landscape and enhancing its utility as a reference.


The following genome assembly pipeline is summarized as below:


![workflow](https://github.com/user-attachments/assets/d1434268-683b-4299-aed8-d8459e65046e)




Souce of data:

PacBioHiFi reads and Hi-C long reads can be retreived from Bioproject ID PRJNA990781. 

Objective

**1. Using hifiasm to perform PacBio genomes aseembly and purge duplicates sequences:**


script to execute at HPC:

`module load miniconda`

`source activate hifiasm`

`hifiasm -o input.asm -t 16 input.fastq`

`awk '/^S/{print ">"$2;print $3}' input.gfa > output.fa`

**2. Use purge_dups to purge haplotigs and overlaps in an assembly based on read depth** 

Purge_dups Installation:

`git clone https://github.com/dfguan/purge_dups.git`
`cd purge_dups/src && make`

Step 1. Use pd_config.py to generate a configuration file.

`nano pb.fofn`

_add abosolute path to fastq_

`~/pd_config.py  -n config.Pregina.json ~/input.fa pb.fofn`

Step 2. install python 3 at HPC

`module load miniconda`

`conda create -n python3 python=3.0`

Step 3. execute the script at HPC

`module load miniconda`

`source activate python3`

`~/run_purge_dups.py -p bash ~/config.input.json ~/purge_dups/src Blowfly`


**3. Using Blobtools to filter the contaminated sequences**

The pipeline insruction and guide can be found in [Blobtool 2](https://blobtoolkit.genomehubs.org/blobtools2/). 

**4. Apply Yahs workflow to perform Hi-C integrated assembly**

*** At first step, use Arima mapping to create bam file of Hi-C fasq to genome mapping file***


**BWA index**



`SRA='HHY-1_'
LABEL='overall_exp_name'
module load bwa
module load samtools
module load jdk
IN_DIR='~/fastq/'
REF='~/input.fa'
FAIDX='$REF.fai'
PREFIX='/~/bwainput/'
RAW_DIR='/rawdata/mapping'
FILT_DIR='/~/filtered/'
FILTER='~/filter_five_end.pl'
COMBINER='~/two_read_bam_combiner.pl'
STATS='~/mapping_pipeline-master/get_stats.pl'
PICARD='/home/slin023/picard/build/libs/picard.jar'
TMP_DIR='~/mapping/tmp_files'
PAIR_DIR='~/mapping/bams'
REP_DIR='~/deduplicated/'
REP_LABEL=$LABEL\_rep1
MERGE_DIR='~/replicates'
MAPQ_FILTER=10
CPU=12`

`echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR`

`echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files`
`bwa index -a bwtsw $REF`

`echo "### Step 1.A: FASTQ to BAM (1st)`
`bwa mem -t $CPU $REF $IN_DIR/HHY-1_R1.fastq | samtools view -@ 12  -o $RAW_DIR/HHY_1_R1.bam`

`echo "### Step 1.B: FASTQ to BAM (2nd)`
`bwa  mem -t $CPU $REF $IN_DIR/HHY-1_R2.fastq | samtools view -@ 12 -o $RAW_DIR/HHY_1_R2.bam`

`echo "### Step 2.A: Filter 5' end (1st)`
`samtools view -h $RAW_DIR/HHY_1_R1.bam | perl $FILTER | samtools view -@ 12 -o $FILT_DIR/HHY_1_R1.bam`

`echo "### Step 2.B: Filter 5' end (2nd)`
`samtools view -h $RAW_DIR/HHY_1_R2.bam | perl $FILTER | samtools view -@ 12  -o $FILT_DIR/HHY_1_R2.bam`

`echo "### Step 3A: Pair reads & mapping quality filter`
`perl $COMBINER $FILT_DIR/HHY_1_R1.bam $FILT_DIR/HHY_1_R2.bam $SAMTOOLS $MAPQ_FILTER | samtools view -bS -t $FAIDX - | samtools sort -@ 12 -o $TMP_DIR/$SRA.bam`

**Step 4A done seperately because combine perl script is looking for abosolute path of samtools**

`perl /~/two_read_bam_combiner.pl ~/filtered/HHY_1_R1.bam ~/mapping/filtered/HHY_1_R2.bam ~/spack/applications/gcc-8.2.0/samtools-1.9-o53igvdgqlbqwbewcsot53hcgwwsvel4/bin/samtools 10 > ~/mapping_pipeline-master/all.bam` 

`samtools view ~/all.bam -t ~/bwa/input.fa.fai  -o ~/mapping/all1.bam`


`samtools sort -@ 12 -n -o $TMP_DIR/$SRA.bam /scratch/mdegenna/slin023/mapping/Arima_mapping.bam`



**5. Run yahs to use Arima mapping .bam as input to create contact matrices **

`export PATH=$PATH:~/yahs`


`yahs ~/input.fa ~/Arima_mapping.bam`

This step generates 3 output files:

1. yahs.out.bin
2. yahs.out_scaffolds_final.agp
3. yahs.out_scaffolds_final.fa

Generate HiC contact map that can be loaded by Juicebox

`export PATH=$PATH:~/yahs`

`juicer pre -a -o out_JBAT ~/all1.bed ~/yahs.out_scaffolds_final.agp  >out_JBAT.log ~/input.fa.fai 2>&1`


**make sure to create fai file for original genome before run the command above**


The step generated 5 output files:

* out_JBAT.liftover.agp
* out_JBAT.assembly.agp
* out_JBAT.assembly
* out_JBAT.txt
* out_JBAT.log

out_JBAT.log will tell you what command you should run next

You need to run:

`module load jdk `

`java -Xmx36G -jar ~/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic <(echo "assembly 534549832")`

This will generate "out_JBAT.hic"

Download Juicer Box locally, use the .agp file generated from step 1 and .hic matrix file to manually edit the genome


If you have output file, use this command to generate the fasta file based on your matrix index

`juicer post -o out_JBAT revised.out_JBAT.assembly ~/yahs/out_JBAT.liftover.agp ~/input.fa`








Reference:
1. Cheng H, Concepcion GT, Feng X, Zhang H, Li H. Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nat Methods. 2021;18: 170–175.
2. Guan D, McCarthy SA, Wood J, Howe K, Wang Y, Durbin R. Identifying and removing haplotypic duplication in primary genome assemblies. Bioinformatics. 2020;36: 2896–2898.
3. Challis R, Richards E, Rajan J, Cochrane G, Blaxter M. BlobToolKit – Interactive Quality Assessment of Genome Assemblies. G3 Genes|Genomes|Genetics. 2020;10: 1361–1374.
4. Zhou C, McCarthy SA, Durbin R. YaHS: yet another Hi-C scaffolding tool. Bioinformatics. Edited by C. Alkan; 2023.
5. Durand NC, Robinson JT, Shamim MS, Machol I, Mesirov JP, Lander ES, et al. Juicebox Provides a Visualization System for Hi-C Contact Maps with Unlimited Zoom. Cell Syst. 2016;3: 99–101.
6. Palmer JM. Funannotate: pipeline for genome annotation. 2016.

