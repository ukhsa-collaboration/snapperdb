__author__ = 'flashton'

'''
This section of the code is heavily influenced by the design of Aaron Quinlan and Nick Loman's poretools package - check it out!
https://github.com/arq5x/poretools

Feature wish list

Want an ipython notebook that speaks to the gdw-sequencing table, you can give it a SNP address and it shows you the
distribution of that snp address over time (and geog).

Also, should work the other way around and given an ebg, tell you the most frequent SNP addresses that week, last week etc.

'''

import argparse
from __init__ import __version__, parse_config
from gbru_vcf import fastq_to_vcf
from snpdb import vcf_to_db, make_snpdb, get_the_snps, update_distance_matrix, qsub_to_check_matrix



class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr", action="store_true", dest="quiet")

def run_command(args):
    config_dict = parse_config(args)

    if args.command == 'fastq_to_db':
        args.fastqs = sorted(args.fastqs)
        vcf = fastq_to_vcf(args, config_dict)
        vcf_to_db(args, config_dict, vcf)

    elif args.command == 'fastq_to_vcf':
        args.fastqs = sorted(args.fastqs)
        fastq_to_vcf(args, config_dict)

    elif args.command == 'vcf_to_db':
        ## third argument is for an instance of a vcf class, which doesnt exist in this case
        vcf_to_db(args, config_dict, None)

    elif args.command == 'update_distance_matrix':
        update_distance_matrix(config_dict, args)

    elif args.command == 'qsub_to_check_matrix':
        qsub_to_check_matrix(config_dict, args)

    elif args.command == 'make_snpdb':
        make_snpdb(config_dict)

    elif args.command == 'get_the_snps':
        print str(args.out)
        get_the_snps(args, config_dict)



def main():
    parser = argparse.ArgumentParser(prog='snapperdb', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed snapperdb version", action="version",
                        version="%(prog)s " + str(__version__))

    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    ## metavar is a variable name for use in help
    ## nargs = '+' means that fastqs will take mutliple files and return a list and return an error if not at least 1

    parser_fastq_to_db = subparsers.add_parser('fastq_to_db', help='Takes fastqs and a config file, produces vcf and then '
                                                                   'inserts variants into the snpdb specified in the config')
    parser_fastq_to_db.add_argument('fastqs', metavar='FASTQ file(s)', nargs='+', help='At least one fastq file')
    parser_fastq_to_db.add_argument('-c', dest='config_file', metavar='Config file', required=True,
                                    help='The name of a config file in the user_configs directory (not the full path)')

    parser_fastq_to_vcf = subparsers.add_parser('fastq_to_vcf', help='Takes fastqs and a config file and produces a vcf and '
                                                                     'serialised SNPs and ignored positions')
    parser_fastq_to_vcf.add_argument('fastqs', metavar='FASTQ file(s)', nargs='+', help='At least one fastq file')
    parser_fastq_to_vcf.add_argument('-c', dest='config_file', metavar='Config file', required=True, help='The name of a config '
                                                                        'file in the user_configs directory (not the full path)')

    parser_vcf_to_db = subparsers.add_parser('vcf_to_db', help='Takes a vcf and a config file, parses the vcf and then adds to '
                                                               'snpdb specified in the config file')
    parser_vcf_to_db.add_argument('vcf', metavar='VCF file', nargs='+', help='A vcf file (generated using '
                                                                                'emit_all_positions?)')
    parser_vcf_to_db.add_argument('-c', dest='config_file', metavar='Config file', required=True, help='The name of a config '
                                                                    'file in the user_configs directory (not the full path)')

    parser_make_snpdb = subparsers.add_parser('make_snpdb', help='Takes a config and makes a snpdb')
    parser_make_snpdb.add_argument('-c', dest='config_file', metavar='Config file', required=True, help='The name of a config '
                                                                    'file in the user_configs directory (not the full path)')
    parser_update_distance_matrix = subparsers.add_parser('update_distance_matrix',
                                                          help='Takes a config and updates the distance matrix in the specified '
                                                               'snpdb')
    parser_update_distance_matrix.add_argument('-c', dest='config_file', metavar='Config file', required=True,
                                               help='The name of a config file in the user_configs directory '
                                                    '(not the full path)')
    parser_update_distance_matrix.add_argument('-m', dest='hpc', default='N', help='This is a PHE only function <int>/N, '
                                                                                   'where int is the number of comparisons you '
                                                                                   'want to do on each core')

    parser_qsub_to_check_matrix = subparsers.add_parser('qsub_to_check_matrix', help='This is only for internal use by snapperdb'
                                                                               ' when update matrix is being run in hpc mode.')
    parser_qsub_to_check_matrix.add_argument('-c', dest='config_file', metavar='Config file', required=True, help='The name of a'
                                                                ' config file in the user_configs directory (not the full path)')
    parser_qsub_to_check_matrix.add_argument('-l', dest='strain_list', required=True, help='The list of all the strains '
                                                                                            ' in the SNPdb')
    parser_qsub_to_check_matrix.add_argument('-s', dest='short_strain_list', required=True, help='The list of all the strains '
                                                                                            'already in the distance matrix')
    parser_qsub_to_check_matrix.add_argument('-u', dest='update_list', required=True, help='The list of all the strains to be '
                                                                                        'added to the distance matrix')

    parser_get_the_snps = subparsers.add_parser('get_the_snps', help='Takes a config file, a list, and a bunch of other flags '
                                                                    'and provides you with snps and more')
    parser_get_the_snps.add_argument('-c', dest='config_file', metavar='Config file', required=True, help='The name of a config '
                                                                        'file in the user_configs directory (not the full path)')
    parser_get_the_snps.add_argument('-l', dest='strain_list', required=True)
    parser_get_the_snps.add_argument('-m', dest='snp_co', required=True, help='SNP cut off (has to be integer), strains more '
                                                                             'than this number of '
                                                                             'SNPs from the reference will be excluded from '
                                                                             'the analysis. A sensible starting point is 3000')

    parser_get_the_snps.add_argument('-o', dest='out', help='Prefix for output, will default to '
                                                           '<date>.<snpdb_name>_vs_<list_name>')
    parser_get_the_snps.add_argument('-a', dest='alignment_type', help='Alignment type (W=whole consensus, A=Accessory, '
                                                                      'C=Core), default is Core', default='C')
    parser_get_the_snps.add_argument('-x', dest='mat_flag', help='Would you like a pairwise distance matrix? Y/N (default = '
                                                                'N)', default='N')
    parser_get_the_snps.add_argument('-v', dest='var_flag', help='Would you like a more detailed list of the variant '
                                                                'attributes? Y/N (default = N)')
    parser_get_the_snps.add_argument('-r', dest='ref_flag', help='Would you like the reference genome in the alignment? Y/N ('
                                                                'default = Y)', default='Y')
    parser_get_the_snps.add_argument('-n', dest='rec_file', help='File with list of positions to ignore due to expected '
                                                                 'recombination', default='N')
    parser_get_the_snps.add_argument('-b', dest='back_flag', help='Would you like a background cluster level? SNP cluster level '
                                                                 'from which to take one representative/N', default='N')
    parser_get_the_snps.add_argument('-e', dest='meta_flag', help='some value from the metadata in strain_stats, '
                                                                  'every strain with this meta-data will be included. '
                                                                  'e.g. (e.g. stx:2a,pt:8,row:value)', default='N')

    parser_update_clusters = subparsers.add_parser('update_clusters', help='Given a config file, updates the SNP clustering '
                                                                           'associated with the SNPdb specified in the config.')
    parser_update_clusters.add_argument('-c', dest='config_file', metavar='Config file', required=True, help='The name of a config '
                                                                        'file in the user_configs directory (not the full path)')

    args = parser.parse_args()
    run_command(args)


if __name__ == "__main__":
    main()






