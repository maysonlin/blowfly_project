# Blowfly_project

# Blowfly Genome Assembly

The script was written by Sheng-Hao (Mayson). The pipeline shown is to generate a contiguous and comprehensive chromosomal scale level of Hi-C intergated genome assembly, providing valuable insights into the genomic landscape and enhancing its utility as a reference.


The following pipeline is summarized as below:


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

The pipeline insruction and guide can be found in [Blobtool 2 https://blobtoolkit.genomehubs.org/blobtools2/]. 
