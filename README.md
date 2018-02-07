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

It is available by pip
```bash
$ pip install snapperdb
```

also from conda
```bash
$ conda install -c tdallman snapperdb
```


and from source:
```bash
$ git clone https://github.com/phe-bioinformatics/snapperdb.git
```


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
however it is best practice to set up the *GASTROSNAPPER_CONFPATH* environment variable to point to your config directory. 

The reference genome should have the .fa suffix and the path defaults to
```
/$PATH/snapperdb/reference_genomes
```
however it is best practice to set up the *GASTROSNAPPER_REFPATH* environment variable to point to your directory of reference genomes.


If the reference is a *de novo* assembly place the original FASTQ files in the reference_genomes directory and it will use them to mask ambiguous mapping regions.  They need to be of the format $reference_genome.R1.fastq.gz $reference_genome.R2.fastq.gz.  If you want your variants to be annotated you can add a Genbank file to this directory with the naming convention reference_genome.gb

We recommend that you *do not* user the default postgres username and password and set up a new user for all your Snapper needs.  This user will need to be a Superuser to create the database as SnapperDB relies on the Postgres extension INTARRAY.

To create the database run the command:

```sh
run_snapperdb.py make_snpdb -c $myconfigname
```

This command will execute the SQL to create the postgres database.  Then it will either simulate reads or use the user supplied reads to map against the reference genome.  Regions of ambiguous mapping (and any variants!) are stored in SNP database.


### Populating a Database with FASTQs

Once a SnapperDB instance has been created you will want to populate it with FASTQs.

SnapperDB can import samples in one by one basis using the **fastq_to_db** command

```sh
run_snapperdb.py fastq_to_db -c $myconfigname $FASTQ1 $FASTQ2
```

If you want to batch load a set of samples.  It is recommended you run **fastq_to_vcf** in parallel followed by **fastq_to_db** in serial.


### Populating a Database from a JSON file

If you have received a JSON export of variants from another SNAPPERDB you can import that sample into your own local database.

```sh
run_snapperdb.py import_json -j $myjsonfile
```

There are two modes of import that can be selected by the **-w** flag. R will simply read from SnapperDB database and return the best match(es) for SNP address. W will write to SnapperDB equivalent to if you were importing a VCF.  R mode is still under development

### Updating the distance matrix

After uploading samples to your SnapperDB instance you can populate the distance matrix with pairwise SNP differences.

```sh
run_snapperdb.py update_distance_matrix -c $myconfigname
```

If you have access to a HPC cluster the **-m** flag can be used to partition your job into.  This parameter expects an integer that will be the number of strains you want update in each qsub job.

### Generating SNP Addresses

Using hierarchical single linkage clustering of the pairwise SNP distances we are able to derive an isolate level nomenclature for each genome sequence.  This allows efficient searching of the population studied as well as automated cluster detection.  By default the SNP address performs single linkage clustering at seven SNP thresholds; 250, 100, 50, 25, 10, 5, 0.

```sh
run_snapperdb.py update_clusters -c $myconfigname
```

Isolates that have artificially low genetic distance to others may led to cluster merges with single linkage clustering.  Artificially low genetic distances might be caused by low quality samples or by samples that are mixed with another isolate.  To combat this problem, for each cluster a Z-Score is calculated for each isolate to ascertain how far that isolates average SNP distance deviates from the mean average SNP distance in that cluster.  If an isolates Z-Score is less that -1.75 it is deemed an outlier and will have to be manually accepted or ignored.

Here is an example of the output of update_clusters where an outlier has been identified. Strain 13816_H14212028301-1 has tripped the Z-Score threshold in the 50 SNP cluster '1', its average distance to other isolates in this cluster was 42.6 which resulted in a Z-Score just less than -1.75.


```
### Cluster level 250 :23:41:32.847220
### Cluster level 100 :23:53:40.710755
### Cluster level 50 :00:00:35.531624
### Cluster level 25 :00:03:01.664428
### Cluster level 10 :00:03:51.182371
### Cluster level 5 :00:04:37.412988
### Cluster level 0 :00:05:26.974864
### Getting previously checked outliers:00:06:26.419278
### Checking Clusters:00:06:26.540735
33-19-8101311 [1, 1, 1, 124, 3121, 4294, 5529]
### Investigating  250 1
### Average Distance  417.962973573
### Standard Deviation  216.131996093
### Investigating  100 1
### Average Distance  69.2875007924
### Standard Deviation  22.0576314
### Investigating  50 1
### Average Distance  67.2382283514
### Standard Deviation  14.066770828
### Outlier Z-Score  13816_H14212028301-1 42.6178690009 -1.75024955276
### Investigating  25 124
### Average Distance  60.1811010179
### Standard Deviation  12.4097892225
### Investigating  10 3121
### Average Distance  10.041015625
### Standard Deviation  3.85288497424
### Investigating  5 4294
### Average Distance  0.0
### Investigating  0 5529
### Average Distance  0.0
```


At this stage it is best to manually inspect the SNP alignment for the outlier (see Generating Alignments for Phylogenetic Analysis below), does it have more 'N' bases then expected?

If you are happy to accept the outlier you can run the following

```sh
run_snapperdb.py accept_outlier -c $myconfigname -n $nameofstrain
```

If you think the outlier is best ignored from this and future analysis then run

```sh
run_snapperdb.py ignore_isolate -c $myconfigname -n $nameofstrain
```

Now you can run update_clusters again and generate the SNP addresses

### Generating Alignments for Phylogenetic Analysis

The command **get_the_snps** is used to produce alignments to use for phylogenetic analysis.  In it's simplest form the **get_the_snps** function takes the name of the config file and list of strains in your database that you want to generate the alignment from.

```sh
run_snapperdb.py get_the_snps -c $myconfigname -l $mylistofstrains
```

There are lots of other options in **get_the_snps**.  To see them supply **-h**.  We think we have provided fairly sensible defaults.  One key option is **-m** which sets the maximum number of SNPs away from the reference that's allowed (Default 5000).  Strains with more SNPs than this will be ignored.  Another important option is **-a** which is alignment type.  The three options available are; **C** (Core) -  only allows positions that are present in all strains, **A** (Soft Core) - produces alignments where at least specified percentage of the samples are A/C/T/G at each position (Default A:80), **W** produce whole genome alignments.  If you wish to mask regions of the genome out of the alignment, for example parts of the reference known to have recombined the you can either supply a file of coordinates to ignore with the **-n** flag or supply a GFF file produced by gubbins for example with the **-ng** option.  The **-b** flag allows the list of strains supplied with the -l flag to be complemented with other background strains from the database.  To include a representative from each 100 SNP cluster you would provide the option **t100**.  Finally if you would like to output a list of annotated variants or a distance matrix of pairwise SNP distances provide a 'Y' option to the **-v** and **-x** flags respectively.

### Exporting variants into a JSON object for importing into another database

This will create a set of JSON files containing all the information to share to another database

```sh
run_snapperdb.py export_json -c $myconfigname -l $mylistofstrains
```
