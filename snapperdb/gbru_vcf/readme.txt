
This module will take

* fastq(s)
* a config file (containing reference genome name).

This module will produce

* sorted.bam.gz (in tmp)
* parsed depth file (in tmp)
* vcf (if depth passes)
* vcf.pickle (tmp)

Dependencies for this module include

* bwa
* gatk
