import os
import logging
from vcf import Vcf

__author__ = 'flashton'

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
    logger.info('Making sorted bam')
    fastq_bam_vcf.make_sorted_bam(args)
    logger.info('Making vcf')
    fastq_bam_vcf.make_vcf(args)
    logger.info('Parsing vcf')
    fastq_bam_vcf.read_vcf()
    logger.info('Checking VCF length')
    fastq_bam_vcf.check_len_vcf(config_dict)
    if args.command == 'fastq_to_vcf':
        logger.info('Pickling variants and ignored positions')
        fastq_bam_vcf.pickle_variants_and_ignored_pos()
    elif args.command == 'fastq_to_db':
        logger.info('Pickling variants and ignored positions')
        fastq_bam_vcf.pickle_variants_and_ignored_pos()
        logger.info('Returning instance of Vcf class')
        return fastq_bam_vcf

def fastq_to_vcf_multi_contig(args, config_dict):
    print 'running fastq_to_vcf_multi_contig'
    fastq_bam_vcf = Vcf()
    fastq_bam_vcf.parse_config_dict(config_dict)
    fastq_bam_vcf.define_class_variables_and_make_output_files(args)
    fastq_bam_vcf.make_sorted_bam(args)
    fastq_bam_vcf.make_vcf(args)
    vcf_container = fastq_bam_vcf.read_multi_contig_vcf()
    fastq_bam_vcf.check_len_vcf(config_dict)
    if args.command == 'fastq_to_vcf':
        fastq_bam_vcf.pickle_variants_and_ignored_pos()


def parse_vcf_for_mixed(args, config_dict):
    if args.vcf_file.endswith('.gz'):
        os.system('gunzip {0}'.format(args.vcf_file))
        args.vcf_file = os.path.splitext(args.vcf_file)[0]
    vcf = Vcf()
    vcf.parse_config_dict(config_dict)
    vcf.sample_name = os.path.basename(args.vcf_file[0]).split(os.extsep)[0]
    vcf.ad_cutoff = float(args.ad_ratio)
    vcf.vcf_filehandle = args.vcf_file
    # vcf.read_rec_file(args.rec_file)
    vcf.read_vcf()
    print 'There are ', vcf.number_mixed_positions, ' mixed positions.'
    with open('{0}/{1}.positions.pick'.format(args.outdir, vcf.sample_name), 'wb') as fo:
        pickle.dump(vcf.mixed_positions, fo)
        # pickle.dump(vcf.number_mixed_positions)
    # with open('{0}/{1}.positions.txt'.format(args.outdir, vcf.sample_name), 'w') as fo:
    #    fo.write('{0}\t{1}\n'.format(vcf.sample_name, vcf.number_mixed_positions))
    # os.system('gzip {0}'.format(args.vcf_file))
