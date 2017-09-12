import os
import logging
from vcf import Vcf

__author__ = 'gidis'

try:
    import cPickle as pickle
except:
    import pickle

def fastq_to_vcf(args, config_dict):
    logger = logging.getLogger('snapperdb.fastq_to_vcf')
    logger.info('Running fastq_to_vcf')
    fastq_bam_vcf = Vcf()
    logger.info('Parsing config_dict')
    fastq_bam_vcf.parse_config_dict(config_dict)
    logger.info('Defining class variables and making output files')
    fastq_bam_vcf.define_class_variables_and_make_output_files(args)
    logger.info('Running Pheonix')
    fastq_bam_vcf.run_phoenix(args)
    return fastq_bam_vcf

    
def make_fastq(args, config_dict):

    logger = logging.getLogger('snapperdb.fastq_to_vcf')
    logger.info('Running fastq_to_vcf')
    fastq_bam_vcf = Vcf()
    logger.info('Parsing config_dict')
    fastq_bam_vcf.parse_config_dict(config_dict)
    logger.info('Defining class variables and making output files')
    fastq_bam_vcf.define_fastq_paths(args)
    fastq_bam_vcf.define_class_variables_and_make_output_files(args)
    logger.info('Making FASTQs')
    fastq_bam_vcf.make_ref_fastqs(args)
    logger.info('Running Pheonix')
    fastq_bam_vcf.run_phoenix(args)
    return fastq_bam_vcf
