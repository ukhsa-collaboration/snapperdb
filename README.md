# SnapperDB README

Welcome to SnapperDB, congratulations on having the excellent taste to use our software for all your FASTQ -> phylogeny needs. 

First things first, we are going to assume you work at PHE on the HPC infrastructure. If you don't then I need to update this documentation.

### Version

0.2

### Basic Useage

You don't need to clone the git repo as the release branch of SnapperDB is available as a module on the cluster.

All you need to do is

```sh
$ module load snapperdb/0.2
```

After that, the SnapperDB_main.py executable will be in your path. 
```sh
$ SnapperDB_main.py 
```
SnapperDB obeys Torsten's first law - 'print something when no parameters supplied'. So you should see some info about useage printed following the above command. More info available if you ask for it.

```sh
$ SnapperDB_main.py 
```

### Typical Workflow

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

Once you have made your config file and filled it with your info, save it in as e.g. *gas_config.txt* 
```
/phengs/hpc_software/snapperdb/0.2/user_configs
```
Now you are ready to make your SNPdb. Run 

```sh
SnapperDB_main.py make_snpdb -c gas_config.txt
```
Hopefully, SnapperDB will read your config files and create a database for you on the postgres server you gave it.

Then, save your reference genome (with bwa index, [samtools faidx and picard sequnced dictionary for vcf calling]) in: 
```
/phengs/hpc_software/snapperdb/0.2/reference_genomes
```

Now, you are ready to run some samples!

There are two intended modes of operation, sample-by-sample and high-throughput.

### Sample-by-sample (s-b-s)

If you are not running anything in parallel, you can run s-b-s.

```sh
SnapperDB_main.py fastq_to_db -c gas_config.txt reads.1.fq reads.2.fq
```
This will run the whole pipeline and deposit the variants and ignored positions in the SNPdb specificed in your config file.

### High-thoughput (h-t)

If you want to run a large number of samples in parallel (i.e. via qsub on the cluster) you should use h-t. This essentially divides the workflow into two, fastq_to_vcf (which is slow and can be done in parallel) and vcf_to_db (which is fast and must be done serially).

Launch lots of qsub jobs that run:

```sh
SnapperDB_main.py fastq_to_vcf -c gas_config.txt reads.1.fq reads.2.fq
```

Then, either set up a qsub wait job or, just manually wait for all the fastq_to_vcf jobs to finish and then run a single job that will run:
```sh
SnapperDB_main.py vcf_to_db -c gas_config.txt sample.vcf
```
on each vcf you are interested in (in serial!).

### Get the SNPs!

Now, you should have lots of lovely data in your SNPdb, get it out by

```sh
SnapperDB_main.py get_the_snps -c gas_config.txt -l strain_list
```
where the strain_list is a list of strain ids that are in your SNPdb. There are lots of other options in get_the_snps. To see them, supply -h. There are sensible-ish defaults for most of these. One of the most relevant to your final phylogeny of all of these is ALIGNMENT_TYPE (-a), which defaults to core. This means that only positions that were present in each isolate in your analysis will be included in your alignment. If you set this to A (i.e. -a A) then SNPs in positions that were only present in a single isolate will be included. Core is the default because it is more robust to strains that have a large number of SNPs in spurious positions. However, Accessory is a very useful option as it is more robust to strains that ignore positions that are true SNPs, that contribute to your phylogeny, possibly because the original isolate was a mix of two strains, leading to bad AD ratios for those positions. With core, if they are ignored in one strain, they are ignored in all strains and so your tree will collapse.

get_the_snps will output you a fasta file that is suitable for running through e.g. FasTree or Raxml.


