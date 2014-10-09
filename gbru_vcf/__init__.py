__author__ = 'flashton'

import os
from gbru_vcf_modules import Vcf


def fastq_to_vcf(opts, config_dict):
    fastqs = sorted(opts.fastqs)
    fastq_bam_vcf = Vcf()
    fastq_bam_vcf.sample_name = os.path.basename(fastqs[0]).split(os.extsep)[0]
    fastq_bam_vcf.make_sorted_bam(opts)
    #fastq_bam_vcf.make_vcf(opts)
    #return fastq_bam_vcf
