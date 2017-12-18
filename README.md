# SnapperDB README


Welcome to SnapperDB, A database solution for routine sequencing analysis of bacterial isolates.

SnapperDB is a python application that sits upon one or more postgres databases to manage reference based SNP typing of bacterial isolates.

SnapperDB can take a pair of Illumina FASTQ sequencing reads and execute a user-defined variant calling pipeline storing the resultant variant calls and absent positions for each isolate.

As the database is populated a pair-wise distance matrix of SNP distances is calculated that can be used to generate an isolate level hierarchical clustering nomenclature - the SNP Address.

SnapperDB instances can be queried to produce alignments for phylogenetic analysis.

SnapperDB has been used internally within Public Health England to process >20,000 isolates of *Salmonella*, *E. coli* and other gastrointestinal pathogens.

SnapperDB is intended to be used at the clonal complex / eBURST group level where isolates are less than 10K SNPs from the reference genome.

---

### Version

1.0.1

---

### Installation:

SnapperDB is available at https://github.com/phe-bioinformatics/snapperdb

From source:

```bash
$ git clone https://github.com/phe-bioinformatics/snapperdb.git
$ pip2 install -e snapperdb
```

Add the installation destination to your PYTHONPATH or run from the installation path.

---

### Dependencies:

**Python**

- Python >= 2.7
- biopython
- psycopg2
- paramiko
- hashids
- joblib

**Postgres**

Postgres package can be downloaded from https://www.postgresql.org/

**PHEnix**

PHEnix is available from https://github.com/phe-bioinformatics/PHEnix

**Samtools**

Samtools can be downloaded from https://github.com/samtools/samtools. It is used to filter and convert to SAM/BAM files and in mpileup variant caller.

**Picard**

The Picard tool suite is available from http://broadinstitute.github.io/picard/

Set *PICARD_JAR* - full path to the GATK Java archive.

**GATK**

GATK is available from https://www.broadinstitute.org/gatk/. Please read the licensing information before using (https://www.broadinstitute.org/gatk/about/#licensing)

Set *GATK_JAR* - full path to the GATK Java archive.

**BWA**

The BWA mapper can be downloaded from http://bio-bwa.sourceforge.net/.

---

### Creating a Database

The first step of SnapperDB is to create a database to populate.  To do this a reference genome is required in FASTA format and config file.
The config file is a tab delimited file containing information about the connection, the name of the reference genome your choice of mapper and variant caller and some user defined features.   Currently SnapperDB supports BWA and GATK only.

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
mapper_threads 4
variant_caller gatk
variant_caller_threads 4

```

*depth_cutoff* is the minimum consensus depth, positions below this will be ignored.

*mq_cutoff* is the mapping quality below which a position will be ignored

*ad_cutoff* is the proportion of the majority variant below which a position is ignored

*average_depth_cutoff* is the average coverage across the genome, samples with coverage below this will be ignored.

The default home for the configs is in
```
/$PATH/snapperdb/user_configs
```

The reference genome should have the .fa suffix and be placed in.
```
/$PATH/snapperdb/reference_genomes
```

If the reference is a *de novo* assembly place the original FASTQ files in the reference_genomes directory and it will use them to mask ambiguous mapping regions.  They need to be of the format $reference_genome.R1.fastq.gz $reference_genome.R2.fastq.gz.  If you want your variants to be annotated you can add a Genbank file to this directory with the naming convention reference_genome.gb


To create the database run the command:

```sh
snapperdb.py make_snpdb -c $myconfigname
```

This command will execute the SQL to create the postgres database.  Then it will either simulate reads or use the user supplied reads to map against the reference genome.  Regions of ambiguous mapping (and any variants!) are stored in SNP database.


### Populating a Database with FASTQs

Once a SnapperDB instance has been created you will want to populate it with FASTQs.

SnapperDB can import samples in one by one basis using the **fastq_to_db** command

```sh
snapperdb.py fastq_to_db -c $myconfigname $FASTQ1 $FASTQ2
```

If you want to batch load a set of samples.  It is recommended you run **fastq_to_vcf** in parallel followed by **fastq_to_db** in serial.


### Populating a Database from a JSON file

If you have received a JSON export of variants from another SNAPPERDB you can import that sample into your own local database.

```sh
snapperdb.py import_json -j $myjsonfile
```

There are two modes of import that can be selected by the **-w** flag. R will simply read from SnapperDB database and return the best match(es) for SNP address. W will write to SnapperDB equivalent to if you were importing a VCF.  R mode is still under development

### Updating the distance matrix

After uploading samples to your SnapperDB instance you can populate the distance matrix with pairwise SNP differences.

```sh
snapperdb.py update_distance_matrix -c $myconfigname
```

If you have access to a HPC cluster the **-m** flag can be used to partition your job into.  This parameter expects an integer that will be the number of strains you want update in each qsub job.

### Generating SNP Addresses

Using hierarchical single linkage clustering of the pairwise SNP distances we are able to derive an isolate level nomenclature for each genome sequence.  This allows efficient searching of the population studied as well as automated cluster detection.  By default the SNP address performs single linkage clustering at seven SNP thresholds; 250, 100, 50, 25, 10, 5, 0.

```sh
snapperdb.py update_clusters -c $myconfigname
```

### Generating Alignments for Phylogenetic Analysis

The command **get_the_snps** is used to produce alignments to use for phylogenetic analysis.  In it's simplest form the **get_the_snps** function takes the name of the config file and list of strains in your database that you want to generate the alignment from.

```sh
SnapperDB_main.py get_the_snps -c $myconfigname -l $mylistofstrains
```

There are lots of other options in **get_the_snps**.  To see them supply **-h**.  We think we have provided fairly sensible defaults.  One key option is **-m** which sets the maximum number of SNPs away from the reference that's allowed (Default 5000).  Strains with more SNPs than this will be ignored.  Another important option is **-a** which is alignment type.  The three options available are; **C** (Core) -  only allows positions that are present in all strains, **A** (Soft Core) - produces alignments where at least specified percentage of the samples are A/C/T/G at each position (Default A:80), **W** produce whole genome alignments.  If you wish to mask regions of the genome out of the alignment, for example parts of the reference known to have recombined the you can either supply a file of coordinates to ignore with the **-n** flag or supply a GFF file produced by gubbins for example with the **-ng** option.  The **-b** flag allows the list of strains supplied with the -l flag to be complemented with other background strains from the database.  To include a representative from each 100 SNP cluster you would provide the option **t100**.  Finally if you would like to output a list of annotated variants or a distance matrix of pairwise SNP distances provide a 'Y' option to the **-v** and **-x** flags respectively.

### Exporting variants into a JSON object for importing into another database

This will create a set of JSON files containing all the information to share to another database

```sh
SnapperDB_main.py export_json -c $myconfigname -l $mylistofstrains
```
