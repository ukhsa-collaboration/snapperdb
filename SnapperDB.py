__author__ = 'flashton'

import argparse

from gbru_vcf import fastq_to_vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('command', help='This can be either:\tfastq_to_vcf\tfastq_to_db\tvcf_to_db\n')
    parser.add_argument('config', help='You always need to pass a config')


    ftv_group = parser.add_argument_group('fastq_to_vcf/fastq_to_db', 'In addition to a config, required arguments for '
                                                                      'fastq_to_vcf or fastq_to_db are:')
    #ftv_group.add_argument('--config')
    ftv_group.add_argument('--fastq_1')
    ftv_group.add_argument('--fastq_2')

    vtd_group = parser.add_argument_group('vcf_to_db', 'In addition to a config, required arguments for vcf_to_db are:')
    #vtd_group.add_argument('--config')
    vtd_group.add_argument('--vcf')


    opts = parser.parse_args()

    print opts.command
    print opts.config
    print opts.fastq_1
    print opts.fastq_2

    if opts.command == 'fastq_to_db':


        pass

