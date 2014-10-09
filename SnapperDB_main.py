__author__ = 'flashton'

'''
Code heavily influenced by Aaron Quinlan and Nick Loman's poretools package - check it out!
https://github.com/arq5x/poretools
'''

import argparse
import os
from __init__ import __version__, parse_config
from gbru_vcf import fastq_to_vcf



class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr", action="store_true", dest="quiet")

def run_command(args):
    config_dict = parse_config(args)
    if args.command == 'fastq_to_db':
        fastq_to_vcf(args, config_dict)

    elif args.command == 'fastq_to_vcf':
        pass
    elif args.command == 'vcf_to_db':
        pass

def main():
    parser = argparse.ArgumentParser(prog='SnapperDB', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed SnapperDB version",
                        action="version",
                        version="%(prog)s " + str(__version__))

    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    ## metavar is a variable name for use in help
    ## nargs = '+' means that fastqs will take mutliple files and return a list and return an error if not at least 1

    parser_fastq_to_db = subparsers.add_parser('fastq_to_db', help='Takes fastqs and a config file, produces vcf and then '
                                                                   'inserts variants into the snpdb specified in the config')
    parser_fastq_to_db.add_argument('fastqs', metavar='FASTQ file(s)', nargs='+', help='At least one fastq file')
    parser_fastq_to_db.add_argument('-c', dest='config_file', metavar='Config file', nargs='+', help='A config file in the '
                                                                                             'user_configs directory')

    parser_fastq_to_vcf = subparsers.add_parser('fastq_to_vcf', help='Takes fastqs and a config file and produces a vcf and '
                                                                     'serialised SNPs and ignored positions')
    parser_fastq_to_vcf.add_argument('fastqs', metavar='FASTQ file(s)', nargs='+', help='At least one fastq file')
    parser_fastq_to_vcf.add_argument('-c', dest='config_file', metavar='Config file', nargs='+', help='A config file in the '
                                                                                             'user_configs directory')
    parser_vcf_to_db = subparsers.add_parser('vcf_to_db', help='Takes a vcf and a config file, parses the vcf and then adds to '
                                                               'snpdb specified in the config file')
    parser_vcf_to_db.add_argument('vcf', metavar='VCF file', nargs='+', help='A vcf file (generated using '
                                                                                'emit_all_positions?)')
    parser_vcf_to_db.add_argument('-c', dest='config_file', metavar='Config file', nargs='+', help='A config file in the '
                                                                                             'user_configs directory')

    args = parser.parse_args()

    run_command(args)


if __name__ == "__main__":
    main()