# SnapperDB README


Welcome to SnapperDB, A scalable database for routine sequencing of bacterial isolates.

SnapperDB is a python application that sits upon one or more postgres databases to manage reference based SNP typing of bacterial isolates.

SnapperDB can take a pair of Illumina FASTQ sequencing reads and execute a user-defined variant calling pipeline storing the result variant calls and absent positions for each isolate.

As the database is populated a pair-wise distance matrix of SNP distances is calculated that can be used to generate a isolate level hierichical clustering nomenclature - the SNP Address.

SnapperDB has been used internally within Public Health England to process >20,000 isolates of *Salmonella*, *E. coli* and other gastrointestinal pathogens.

---

### Version

1.0

---

### Installation:

SnapperDB is available at https://github.com/phe-bioinformatics/snapperdb


---

### Dependencies:

**Python**

- Python >= 2.7
- biopython
- pyscopg2

**Postgres**

Postgres package can be downloaded from https://www.postgresql.org/

**PHEnix**

PHEnix is available from https://github.com/phe-bioinformatics/PHEnix

**Samtools**

Samtools can be downloaded from https://github.com/samtools/samtools. It is used to filter and convert to SAM/BAM files and in mpileup variant caller.

**Picard**

The Picard tool suite is available from http://broadinstitute.github.io/picard/

**GATK**

GATK is available from https://www.broadinstitute.org/gatk/. Please read the licencing information before using (https://www.broadinstitute.org/gatk/about/#licensing)

Set *GATK_JAR* - full path to the GATK Java archive.

**BWA**

The BWA mapper can be downloaded from http://bio-bwa.sourceforge.net/.

---

### Creating a Database

The first step of SnapperDB is to create a database to populate.  To do this a reference genome is required in FASTA format and config file.
The config file is a tab deliminated file containing information about the connection, the name of the reference genome your choice of mapper and variant caller and some user defined features.   Currently SnapperDB supports BWA and GATK only.

An example of *Salmonella* Enteritidis is shown below - lets call this file ebg4_config.txt

```
snpdb_name ebg_4_snps
reference_genome AM933172
pg_uname timdallman
pg_pword
pg_host localhost
depth_cutoff 10
mq_cutoff 30
ad_cutoff 0.9
average_depth_cutoff 30
mapper bwa
variant_caller gatk

```

*depth_cutoff* is the minimum consensus depth, positions below this will be ignored.

*mq_cutoff* is the mapping quality below which a position will be ignored

*ad_cutoff* is the proportion of the majority variant below which a position is ignored

*average_depth_cutoff* is the average coverage accross the genome, samples with coverage below this will be ignored.

The default home for the configs is in 
```
/$PATH/snapperdb/user_configs
```

The reference genome should have the .fa suffix and be placed in.
```
/$PATH/snapperdb/reference_genomes
```

If the reference is a *de novo* assembly place the orginal FASTQ files in the reference_genomes directory and it will use them to mask ambigous mapping regions.  They need to be of the format $reference_genome.R1.fastq.gz $reference_genome.R2.fastq.gz.  If you want your variants to be annotated you can add a Genbank file to this directory with the naming convention reference_genome.gb


To create the database run the command:

```sh
snapperdb.py make_snpdb -c $myconfigname
```

This command will execute the SQL to create the postgres database.  Then it will either simulate reads or use the user supplied reads to map against the reference genome.  Regions of ambigous mapping (and any variants!) are stored in SNP database. 


### Populating a Database with FASTQs

Once a SnapperDB instance has been created you will want to populate it with FASTQs.

SnapperDB can import samples in one by one basis using the **fastq_to_db** command

```sh
snapperdb.py fastq_to_db -c $myconfigname $FASTQ1 $FASTQ2
```

If you want to batch load a set of samples.  It is recommended you run **fastq_to_vcf** in parralel followed by **fastq_to_db** in serial. 


### Populating a Database with JSON Export

*To follow*


### Updating the distance matrix

After uploading samples to your SnapperDB instance you can populate the distance matrix with pairwise SNP differnces.

```sh
snapperdb.py update_distance_matrix -c $myconfigname
```

### Generating SNP Addresses

Using hierarchical single linkage clustering of the pairwise SNP distances we are able to derive an isolate level nomenclature for each genome sequence.  This allows efficient searching of the population studied as well as automated cluster detection.  By default the SNP address performs single linkage clustering at seven SNP thresholds; 250, 100, 50, 25, 10, 5, 0.

```sh
snapperdb.py update_clusters -c $myconfigname
```

### Generating Alignments for Phylogenetic Analysis

The command **get_the_snps** is used to produce alignments  





(recombination)

So, you have successfully loaded the module and run SnapperDB_main.py. What next?

First, you need to make a SNPdb to fill with all your lovely data.

Read the output of
```sh
$ SnapperDB_main.py make_snpdb -h
```

Apparently you need a config file. What do they look like? Well, like this:

```
snpdb_name what_you_want_the_snpdb_to_be_called
reference_genome the_reference_genome_name
pg_uname your_postgres_username
pg_pword your_postgres_password
pg_host your_postgres_server
depth_cutoff 10
mq_cutoff 30
ad_cutoff 0.9
```
*depth_cutoff* is the minimum consensus depth, positions below this will be ignored.
*mq_cutoff* is the mapping quality below which a position will be ignored
*ad_cutoff* is the proportion of the majority variant below which a position is ignored


### Get the SNPs!

Now, you should have lots of lovely data in your SNPdb, get it out by

```sh
SnapperDB_main.py get_the_snps -c gas_config.txt -l strain_list
```
where the strain_list is a list of strain ids that are in your SNPdb. There are lots of other options in get_the_snps. To see them, supply -h. There are sensible-ish defaults for most of these. One of the most relevant to your final phylogeny of all of these is ALIGNMENT_TYPE (-a), which defaults to core. This means that only positions that were present in each isolate in your analysis will be included in your alignment. If you set this to A (i.e. -a A) then SNPs in positions that were only present in a single isolate will be included. Core is the default because it is more robust to strains that have a large number of SNPs in spurious positions. However, Accessory is a very useful option as it is more robust to strains that ignore positions that are true SNPs, that contribute to your phylogeny, possibly because the original isolate was a mix of two strains, leading to bad AD ratios for those positions. With core, if they are ignored in one strain, they are ignored in all strains and so your tree will collapse.

get_the_snps will output you a fasta file that is suitable for running through e.g. FasTree or Raxml.


