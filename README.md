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

1.0.6

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

### Software Dependencies and Environment Variables:

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

Implementation does not supoprt GATK v4 as UnifiedGenotyper has been depreciated.

Set *GATK_JAR* - full path to the GATK Java archive.

**BWA**

The BWA mapper can be downloaded from http://bio-bwa.sourceforge.net/.

---

### Functionality 

SnapperDB has several sub commands that can be invoked and are summarised below.

| Command | Purpose |
|---------|---------|
|make_snpdb | Create a new SnapperDB database
|fastq_to_db | Perform reference mapping for a pair of FASTQs, filter the VCF and add to the database
|fastq_to_vcf | Perform reference mapping for a pair of FASTQs and filter VCF
|vcf_to_db | Parse filtered VCF and add to database
|update_distance_matrix | Perform pairwise SNP distance calculations
|update_clusters | Perform single linkage clustering on SNP distance matrix
|ignore_isolate | Remove an isolate from the database
|accept_outlier | Allow an isolate to persist in the database that has failed the outlier test
|get_the_snps | Generate a sequence alignment for phylogenetic inference
|get_strains | Return the strains in the database and their SNP addresses
|export_json | Produce JSON for import to another SnapperDB database
|import_json | Import JSON exported from another SnapperDB database
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


*snpdb_name* is the name of the postgres database

*reference_genome* is the name of the reference isolate (dont include the .fa or .fasta extension)

*pg_uname* is the postgres username (to create a database the user needs to be a postgres SUPERUSER)

*pg_pword* is the postgres password for this user

*pg_host* is the postgres host (set to localhost if the postgres server is running on your local machine)

*depth_cutoff* is the minimum consensus depth, positions below this will be ignored.

*mq_cutoff* is the mapping quality below which a position will be ignored

*ad_cutoff* is the proportion of the majority variant below which a position is ignored

*average_depth_cutoff* is the average coverage across the genome, samples with coverage below this will be ignored.

*mapper* is the mapping software to deploy

*mapper_threads* is the number of threads to use for mapping

*variant_caller* is the variant calling software to deploy

*variant_caller_threads* is the number of threads to use for variant calling


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


If the reference is a *de novo* assembly place the original FASTQ files in the reference_genomes directory and it will use them to mask ambiguous mapping regions.  They need to be of the format $reference_genome.R1.fastq.gz $reference_genome.R2.fastq.gz.  If you want your variants to be annotated you can add a Genbank file to this directory with the naming convention reference_genome.gbk

If you want to see examples of the Public Health England SnapperDB configs and reference genomes please visit https://github.com/phe-bioinformatics/snapperdb_references


We recommend that you *do not* use the default postgres username and password and instead set up a new user for all your SnapperDB needs.  This user will need to be a postgres Superuser to create the database as SnapperDB relies on the Postgres extension INTARRAY.  If you need to grant your SnapperDB user superuser privelages this is best done through the PSQL interface with the command:

```sh
ALTER USER pg_uname WITH SUPERUSER
```


Subsequent users can be set up to read and write to a SnapperDB instance which does not need superuser privileges with the PSQL command:

```sh
GRANT ALL on database snpdb_name to pg_uname
```

Now everything is set up to create the database execute the run_snapperdb.py python program:

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


This command will create a folder called **snpdb** which will contain the BAM file and the filtered VCF file.

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

The command **get_the_snps** is used to produce alignments to use for phylogenetic analysis.  In it's simplest form the **get_the_snps** function takes the name of the config file and a file containing a list of strains in your database that you want to generate the alignment from, each seperated by a new line.

```sh
run_snapperdb.py get_the_snps -c $myconfigname -l $mylistofstrains
```

There are lots of other options in **get_the_snps**.  To see them supply **-h**.  We think we have provided fairly sensible defaults.  One key option is **-m** which sets the maximum number of SNPs away from the reference that's allowed (Default 5000).  Strains with more SNPs than this will be ignored.  Another important option is **-a** which is alignment type.  The three options available are; **C** (Core) -  only allows positions that are present in all strains, **A** (Soft Core) - produces alignments where at least specified percentage of the samples are A/C/T/G at each position (Default A:80), **W** produce whole genome alignments.  If you wish to mask regions of the genome out of the alignment, for example parts of the reference known to have recombined the you can either supply a file of coordinates to ignore with the **-n** flag or supply a GFF file produced by gubbins for example with the **-ng** option.  The **-b** flag allows the list of strains supplied with the **-l** flag to be complemented with other background strains from the database.  To include a representative from each 100 SNP cluster you would provide the option **t100**.  Finally if you would like to output a list of annotated variants or a distance matrix of pairwise SNP distances provide a 'Y' option to the **-v** and **-x** flags respectively.

### Exporting variants into a JSON object for importing into another database

This will create a set of JSON files containing all the information to share to another database

```sh
run_snapperdb.py export_json -c $myconfigname -l $mylistofstrains
```


### Retrieving strains from a database

If you want to retrieve a list of strains that are in the database and their corresponding SNP Address the command **get_strains** can be used.

```sh
run_snapperdb.py get_strains -c $myconfigname
```


If you want to restrict the list to a certain SNP cluster this can be provided by the **-t** flag.  The format it expects is for example t5:120 which will return all isolates in the 5 SNP cluster numbered 120.

### Deleting or purging your database.

If you want to empty or delete your SnapperDB instance this is best done directly through postgres through the PSQL interface.

```sh
drop database $mydatabase name
```

```sh
delete from strains_snps, strain_clusters, variants, strain_stats, variants, ignored_pos, dist_matrix, cluster_logs
```

---

### Tutorial for setting up SnapperDB Instance

This section details an example of setting up a SnapperDB database for *Salmonella* Enteritidis.

The first thing to do is to select a *Salmonella* Enteritidis reference genome.  At PHE we use the PT4 epidemic strain P125109 - Accession AM933172.  We've have put all the reference genomes we use in the following github repository https://github.com/phe-bioinformatics/snapperdb_references.  You can download the file AM933172.fa from there.  Next we need to create the corresponding config file to tell SnapperDB where to make the postgres database, and how to perform the mapping and variant calling.  The SnapperDB configs can also be found in the above github repo.  The *S* Enteritidis config is called ebg4_config.txt (EBG4 is the MLST Eburst Group that defines *S* Enteritidis).  Download the file and lets look inside:

```sh
snpdb_name ebg_4_snps
reference_genome AM933172
pg_uname postgres
pg_pword postgres
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

snpdb_name is the name of the postgres database we are going to create and we've called it ebg_4_snps.  We've also specified that we want to use the AM933172.fa reference genome.  Note that we don't specify the reference genomes FASTA suffix.  You will also need to supply your Postgres credentials and as this user is going to create a database with some postgres extensions we will need to have superuser privileges.

Before we try and create the database we need to setup some environment variables to tell SnapperDB where to find the references and configs.  To do that on the command line:
```sh
export GASTROSNAPPER_CONFPATH='my path to my configs'
export GASTROSNAPPER_REFPATH='my path to my reference genomes'
```
If you don't do this SnapperDB will look in the reference_genomes directory and the user_configs directory in the cloned github repo.

If you want to annotate the variants in terms of what genes they occur in and their consequence (e.g. synonymous, non-synonymous) we need to also supply a genebank annotation of the reference genome.  The genbank file for AM933172 can be found at NCBI here https://www.ncbi.nlm.nih.gov/nuccore/AM933172.  Download it, save it as AM933172.gbk and put it in your path to your reference genomes.

No we are ready to make our SnapperDB instance for *S* Enteritidis.  Run the following command.

```sh
run_snapperdb.py make_snpdb -c ebg4_config.txt
```

The first thing SnapperDB will do is simulate some FASTQ reads for AM933172.fa with Samtools wgsim and map them against the FASTA file using BWA MEM.  Variant calling will be taken care of by GATK and the VCF filtered to call variants with an MQ > 30, mapped by at least 10 reads in which at least 90% of the bases support the variant call.  Positions that don't meet this criteria will also be recorded as 'ignored positions'.  The SNPs and ignored positions will be injected into the postgres database and the reference genome given the SNP address 1.1.1.1.1.1.1

Now we have a SnapperDB database we can add some other *S* Enteritidis genomes.  Let's download some from the PHE Pathogens *Salmonella* NCBI Bioproject PRJNA248792.  Here is a set of 8 genomes.
- ERR2200244
- SRR5055288
- SRR5194193
- SRR5583186
- SRR5815674
- SRR5850014
- SRR5864444
- SRR6131972 

To upload and perform variant calling on SRR6131972 we can perform the following command:

```sh
run_snapperdb.py fastq_to_db -c ebg4_config.txt SRR6131972_1.fastq.gz SRR6131972_2.fastq.gz
```

If you have access to a cluster you might want to perform the VCF generation in parallel and then upload the VCF one at a time.  Uploading of VCF files has to be done in serial to ensure that the same variants are given the same internal IDs in the database.

```sh
run_snapperdb.py fastq_to_vcf -c ebg4_config.txt SRR6131972_1.fastq.gz SRR6131972_2.fastq.gz
```

```sh
run_snapperdb.py vcf_to_db -c ebg4_config.txt SRR6131972.filtered.vcf
```

We now have a database with 8 isolates mapped to our *S* Enteritdis reference genome AM933172.  Next let's cluster the genomes based on their SNP distances.  SnapperDB employs single linkage hierarchical clustering at 7 SNP distance thresholds, 250, 100, 50, 25, 10, 5 and 0 SNPs.  The resultant cluster membership can then be visualised as a SNP Address.  To do this lets compute a distance matrix of SNP distance that will be stored in the database.

```sh
run_snapperdb.py update_distance_matrix -c ebg4_config.txt
```

Once this is completed it will automatically invoke the clustering algorithm to generate SNP addresses.  You can also run this command separately.

```sh
run_snapperdb.py update_clusters -c  ebg4_config.txt
```

We can then return the SNP Address for the isolates in the database.

```sh
run_snapperdb.py get_strains -c  ebg4_config.txt
```

which will return an output like this:

```
AM933172		1.1.1.1.1.1.1
SRR5194193_1	1.1.1.1.1.1.9
SRR5055288_1	1.1.2.2.2.2.2
SRR5815674_1	1.1.3.3.3.3.3
SRR5864444_1	1.1.3.5.5.5.5
SRR5850014_1	1.1.3.5.5.8.8
SRR6131972_1	1.1.4.4.4.4.4
ERR2200244_1	1.1.7.7.7.7.7
SRR5583186_1	2.2.6.6.6.6.6
```

From this we can see that SRR5583186 is in a different t250 cluster then the rest (it's SNP address starts with a 2 instead of 1) and is therefore at least 250 SNPs distant from the other isolates.  All the other strains are in the same t100 cluster as their SNP addresses start with 1.1., this means that each isolate is within 100 SNPs of at least 1 other isolate in that cluster.  SRR5194193 is in the same t5 cluster as the reference genome as they share the first 6 numbers of the SNP address (1.1.1.1.1.1.), where as SRR5864444 and SRR5850014 are in the same 10 SNP cluster (1.1.3.5.5.) as they share the first 5 numbers.  If we wanted to query the database and just return isolates from a specific SNP cluster we can add an extra command to the **get_strains** call.

```sh
run_snapperdb.py get_strains -c  ebg4_config.txt -t t50:3
```
This just returns:
```
SRR5815674_1	1.1.3.3.3.3.3
SRR5864444_1	1.1.3.5.5.5.5
SRR5850014_1	1.1.3.5.5.8.8
```

Finally, let's generate a SNP alignment that we could use to compute a tree.  To do this we shall use the **get_the_snps** function.  This function expects to be provided with a file with the list of strains we want to use, with each strain on a new line.  Lets create a file called list.txt with each of the 8 strain names in it.  Now we can call **get_the_snps**

```sh
run_snapperdb.py get_the_snps -c ebg4_config.txt -l list.txt -o output
```

This will generate the alignment file output.fa

```
>SRR5055288_1
ACGCACCGTCAGAGACCTCGTTATACAGAAGGAGTGGATAGGCGAGCGTTCTGCGCTCAGATTTTCTGCGGCCGCAACATAACCCCAACCGACCGACCAAAATGGATCCACCCAATGTCGTCTTTCCGCCCTCAAAGGAGCCCTTTCTCCGACGCAGCTGCTCAAACTCTACGATAGATTTGGTCGTTAAGGCACCCCGGACCATGGCTGTACGATCGGACGGGCCGGTCACCAATGTGCTTGAGTTGGGGGTTCATAGCCCAGAGGAGGCAAGGGAGCAAGCNTGGGCATTATCCGGAGGTGTAATGAATCTACCCCCTTAGCCTGGCACGGACTGATCTTCAAATCGATGCGTGGGTGTACCCAGGAATTGTCCCCCGGCCAGACCAACGTGGCCGGTAGGGTACGCCGGTCCTACCTTCAGTCTGTAAATTGATGTATAGATGGGGAATCAACAACCCGGGGCCGGTGGATACACTTGTGTGAGGTCGGCCAACGCGTCCCAGTGCTTTTGGCACCAGGACTCGCCGCCTGGCCTGTAGCGGAGCGTGGCTGCCTGCGGTTATAAATTAAGTGCATTCGGAGAGGTTCTACCTCCAAGGCTAGCCACTACAGGCAGTTGTGGGGTATGGCTAGGGAGCTCTTGAGACGCTAAGCTGCTATCGGG
>SRR5815674_1
ACGCACCGTCAGAGTCCTCGTTATACGGAAGGAGTGGATAGGCGAGCGTACTGCGCTCAGATTTTCTTCAGCTGCAACATAACCCCAACCGACCTACCAACTTGGATCCACCCGATGTCGCCTTTCTGTCCTCAAGGGAGCCCCTTTTCCGATTAAGCTGCTCAGACTCTACGATGGGTTTGGTCGTTAAAGCACCTCGGACCATGGCTGTACGGTCAGACGAGACCGTCACCGACGTGCTTGAGTTGGGGGTTCACAACCCAAAGGATGCAAAGGAGCAGGCATGGGCATCATCCGGAGGTGTAATGAATCTACTCCCTCAGCCTGGCGCGGACTAATCTTCAGATCGATGCGCAAGTATACCCAGGTATAGTCCCCCGGCCAGTCCAACGTGGCCGGTAGGGTACGCCGGTCCTACCTTCAGTCTGTAAGTTGATGTATAGATGGGGAATCAACAACCCGGGGCCGGTGGACACACTTGTGTGAGGTCGGACAACACGCCCCCGTGCTTTTGGCACCAGGACTAGCTACCTGGCCTATAGCGGAGCGTAGCGGCCTGCGGGTATAAAATAAAAGGCATCGGAGAGATTTTACCTCCAAAGTCAGCCATTACAGGCAGCTGGTGGTTGCTGCCAGGGAGCTCTTGAGACACTAAGCTTCTATCGGG
>SRR5850014_1
ACGCACCGTCAGAGTCGNCGTTATACGGAAGGAGTGGATAGGCGAGCGTACTGCGCTTATATTTTCTTCAGCTGCAANATAACCCCAACCGACCTACCAACTTGGATCCACCCGATGTCACCTTTCTGTCCTCAAGGGAGCCCCTTTTCCGATTAGGCTGCTTAGACTCAACGATGGGTTTGGACGTTAAGGCACCTCGGACCATGGCTGTACGGTCAGACGAGACCGTCACCGACGTGCTTGAGTTGGGGGTTCACAACCCAAAGGATGCAAAGGAGCAGGCCTGGGCATCATCCGGATGTGTAATGAATCTACTCCCTCAGCCTGGCGCGGACTGATCTACAGATCGATGCGCGAGTGTACCAAGGAATAGTCCCCCGGCCNGTCCAACGTGGCCGTTAGGGTACGCCGGTCCTACCTTCAGTCTGTAAGTTGATGTATAGATGGGGAATCAACAACCCGGGGCCGGTGGACACACTTGTGTGAGGTCGGCCANCACGCCCCCGCGCTTTTGGCACCAGGATTAGCTACCTGGCCTATAGAGGAGCGTAGCGGCCTGCGGGTATAAAATAAATGGCATCGGAAAGATTCTACCTCCAAAGCCAGCCATTACAGGCAGCTGGTGGTTGCTGCCAGGGAGCTCTTGAGACACTAAGCTTGTATCGGG
>SRR5864444_1
ACGCACCGTCAGAGTCGTCGTTATACGGAAGGAGTGGATAGGCGAGCGTACTGCGCTTATATTTTCTNCAGCTGCAACATAACCCCAACCGACCTACCAACTTGGATCCACCCGATGTCACCTTTCTATCCTCAAGGGAGCCCCTTTTCCGATTAGGCTGCTTAGTCTCAACGATGGGTTTGGACGTTAAGGCACCTCGGACCATGGCTGTACGGTCAGACGAGACCGTCACCGACGTGCTTGAGTTGAGGGTTCACAACCCAAAGGATGCAAAGGAGCAGGCCTGGGCATCATCCGGATGTGTAATGAATCTACTCCCTCAGCCTGGCGCGGACTGATCTACAGATCGATGCGCGAGTGTACCAAGGAATAGTCCCCCGGCCAATCCAACGTGGCCGTTAGGGTACGCCGGTCCTACCTTCAGTCTGTAAGTTGATGTATAGATGGGGAATCAACAACCCGGGGCCGGTGGACACACTTGTGTGAGGTCGGCCAACACGCCCCCGCGCTTTTGGCACCAGGATTAACTACCTGGCCTATAGCGGAGCGTAGCGGCCTGCGGGTATAAAATAAATGGCATCGGAAAGATTCTACCTCCAAAGCCAGCCATTACAGGCAGCTGGTGGTTGCTGCCAGGGAGCTCTTGAGACACTAAGCTTGTATCGGG
>SRR5583186_1
TTTTAATTGTGAGAAACGTACGACCAGAGCAAAACTACGATACAGACTGACTACGTTCAGCCGCCTCGGGAACGTGGCGCAAATTTGAGTTTTCGGACGGAACGAGTCTGTTCAACACTGCCCACTCGCCCGTGCAGTGACTCCCCCCATTCCTCATCAAACCGGATCATATGGCGGGGACGATCGCCGCGATGGCCAGCATTGCATACGACTTGTTAGGCAATCACCATTTCAGTACACCCGACGCAGTGACCTGCGGATAGGGAGTGGCACGATCTTTGAGCTAAAGGGTGCAGGNNNAAAGCGCAGGCTGGTTTTACCTAATTATCGCGTGACGAAAGTTGGAGTAGCGTGCGATCGCGCTCTGAAGCTACACTCGAGTTGGACACGTACATTTGGCGAAGCCCATAAACTACATTCCTGACCCTTGAGCCACCATGGGACAAGAGGGGACGTCGTTTAAAGGTTGCATGCATTTCGCTAGACAACCAACTGCAGCACTTTCATAGGCCGAAGCGTTGAGCCCGTTGCATTGCTCGCGACAGGATACGAGGATGCATGTGCGCGGGACGGGTGCCACTTAGGGGGCACCTTGCTCATGACCGAACGCCCAGAATGCCCAGTAGGGACTATCGCAAGACCGGCAATGTATCGGTACTCGGGTTAG
>ERR2200244_1
ACGCACCGTCAGAGACCTCGTTATACGGAAGGGGTGGATTGGTGAGCGTACTGCTCTCCGATTTTCTGCGGCCACAACATGTCCCCAACCGACCGACCAAAATTGATCCACCCAACGTCGCTTTTCCGCCCTCAAAGGAGCCCCTTCTCCGACTCAGCTGNTCAGACTCTGCAATGTGTTTGGTTGTTAAGGCACTCCAGCCCATGGCTATACGGACAAAAGAGCCCGNCNCCAATGTGCTTGAGTTGGGGGTTCACAGCCCAGAGGAGGCAAGGGAGCAGGCCCGGGCATTATCCACAGGTGTAATGAATCTACNNCCTCAGCCTGGCGTAGACTGCTATTCAGATCGATACGCTAGTGTATCCAGGAATTGTCCCCCGGCCAGACCAACGTGGCCGGTAGGGTACGCCGGTCCTACCTTCAGTCTGCACGTTGATGGATAGATGGGGAATCAACAACTCGGGTCCGGTGGACACACTTGTGTGAGGTTGGCCAACGAGCCCCCGTGCTTTTGGCACCAAGACTCGCTGTCTGGTCTGTAGCGGAGCGTGGGGGCCTGCGGGTATAAAATAAGTTCCATCGGAGAGGTTCTACCTCCAAGGCCAGCCACTACAGGCAGCTGGTGGGTACTGCCAGGGAGCTCTTNAGACACTAAGCTTCGATCGGA
>SRR6131972_1
ACGCGCCGTCAGAGACCTCGTTATACGGAAGGGGTGGATAGGCGAGTGTACTGAGCGCAGATTTTCTGCGGCCGCAACATAACCCCAACCGACTGACTAAAATGGAGTCACCTAGCGTCGCCTTTCCGCTCTCAAAAGAGACTCTTCTCCGACTCAGTTGCTCAGACTCTACGATGGGTTTAGTCGTTAAGGCACCCCGGACCATGGCTGTACGGTCAGACGAGCCCGTCACTAATGTGTTTTGGTTGGGAGTTCACAGCCCAGAGAAGGGGAGGGAGCAGGCCTGGGCATTATCCGCGGGTGTAATGAATCTACTCCCTCAGCCCGGAGCGGACTGATATTCAGGTCGATGCGAGAGTGTACCCAGGAATTGTCACCCGTCCAGATCAACGTGGCCAGTAGGATATGCCGGTCCTCCCTTCAGTATGTAAGTTGATGTATAGATGTGAAATCAACAACTCGGGGCCGATGGACCCACTTGCGTGAGGTCGGCCAACGCGCCCCCGTGCTTTTGGCACCAGGACTCGCTGCCTGTCCTGTAGCGAAGCGTGGGGGCCTGCAGGTATAAAATAAGTGCCATCGGAGAGGTTCTACCTCCAAGGCCAGCCACTACAGGCAGCTGGTGGGTACTGCCAGGGAGCTCTTGTGACACTAAGCTTCGATCGGG
>SRR5194193_1
ACGCACCGTCAGAGACCTCGTTGTACGGAAGGAGTGGATAGGCGAGCGTATAGCGCTCAGATTTTCTGCGGCCGCAATATAACCCCAGCCGACCGACCAAAATGGATCCACCCAACGTCGCCTTTCCGCCTTCAAAGGAGCCCCTTCTCCGACTCAGCTGCTCAGACTCTACGATGGGTTTGGTCATTAAGGCACCCCGGACCATGGCTGTACGGTCAGACGAGCCCGTCACCAATGTGCTTGAGTTGGGGGTTCACAGCCCAGAGGAGACAAGGGAGCAGGCCTGGGCATTATCCGGAGGTGTAATGAATCTACTCCCTCAGCCTGGCGCGGACTGATATTCAGATCGATGCTCGAGTGTACCCAAGAATTGTCCCTCGGCCAGACCAACGTGGCCGGTAGGGTACGCCGGTCCTACCTTCAGTCTGTAAGTTGATGTATAGATGGGGAATCAACAACTCGGGGCCGGTGGACACACTTGTGTGAGGTCGGCCAACGCGCCCCCGTGCTTTTGGCACCAGGACTCGCTGCCCGGCCTGTAGCGGAGCGTGGGGGCCTGCGGGTATAAAATAAGTGCCATCGGAGAAGTTCTACCTCTGAGGCCAGCTACTACAGGCAGCTGGTGAGTACTGCCAGGGAGTTCTTGAGACACTAAGCTTCGATCGGG
```

The default alignment option is called A:80 which means that all variants will be returned as long as that position is an A/T/C/G and not a N in at least 80% of the isolates analysed.  This can be altered with the -m option to be a different threshold or to return the whole genome consensus including invariant positions.

**get_the_snps** can also be used to output a list of the variants as well as a distance matrix.  Let's run it again with the following commands:

```sh
run_snapperdb.py get_the_snps -c ebg4_config.txt -l list.txt -x Y -v Y
```

We now have an output.variants file and an output.matrix file.  If we look in the matrix file we can see the pair-wise SNP distances between each strain.

```
SRR5815674_1	SRR5055288_1	91
SRR6131972_1	SRR5055288_1	92
SRR6131972_1	SRR5815674_1	104
SRR6131972_1	SRR5864444_1	110
SRR6131972_1	SRR5583186_1	503
SRR6131972_1	ERR2200244_1	85
SRR6131972_1	SRR5850014_1	107
SRR5864444_1	SRR5055288_1	98
SRR5864444_1	SRR5815674_1	31
SRR5864444_1	SRR5850014_1	6
SRR5583186_1	SRR5055288_1	571
SRR5583186_1	SRR5815674_1	511
SRR5583186_1	SRR5864444_1	521
SRR5583186_1	SRR5850014_1	516
ERR2200244_1	SRR5055288_1	85
ERR2200244_1	SRR5815674_1	96
ERR2200244_1	SRR5864444_1	103
ERR2200244_1	SRR5583186_1	487
ERR2200244_1	SRR5850014_1	100
SRR5850014_1	SRR5055288_1	95
SRR5850014_1	SRR5815674_1	27
SRR5194193_1	SRR5055288_1	62
SRR5194193_1	SRR5815674_1	75
SRR5194193_1	SRR6131972_1	66
SRR5194193_1	SRR5864444_1	81
SRR5194193_1	SRR5583186_1	479
SRR5194193_1	ERR2200244_1	59
SRR5194193_1	SRR5850014_1	77
```
and the output.variants is the annotated variants from the alignment
```
ID 	CONTIG 							POS 	SNP MUTATION	GENE FUNCTION
147	gi|206707319|emb|AM933172.1|	362833	A	SYNONYMOUS	outer membrane fimbrial usher protein
715	gi|206707319|emb|AM933172.1|	368418	A	A110E		possible transmembrane regulator
148	gi|206707319|emb|AM933172.1|	372650	T	SYNONYMOUS	putative outer membrane efflux lipoprotein
716	gi|206707319|emb|AM933172.1|	389866	G	V132G		probable terminal oxidase subunit II
692	gi|206707319|emb|AM933172.1|	413567	T	G152E		D-alanine:D-alanine ligase A
693	gi|206707319|emb|AM933172.1|	434115	T	G28V		proline-specific permease ProY
149	gi|206707319|emb|AM933172.1|	449946	C	H242P		riboflavin biosynthesis protein RibD
150	gi|206707319|emb|AM933172.1|	468453	C	SYNONYMOUS	hpothetical major facilitator family transport protein
151	gi|206707319|emb|AM933172.1|	469778	G	F52V		putative exported protein
152	gi|206707319|emb|AM933172.1|	470758	C	SYNONYMOUS	putative exported protein
153	gi|206707319|emb|AM933172.1|	471763	C	*136Q		conserved hypothetical protein (pseudogene)
154	gi|206707319|emb|AM933172.1|	473898	T	M19I		cytochrome o ubiquinol oxidase C subunit
```

Finally I'm going to use FastTree to produce a tree of these strains.

```sh
FastTree -nt output.fa
```

![tree](images/snapper.png)

---

### Contact

tim.dallman@phe.gov.uk

