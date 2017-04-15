#!/usr/bin/env python

__author__ = 'gidis'

#imports
import argparse
import logging
import os
import sys
import datetime
from snapperdb import __version__, parse_config
from snapperdb.gbru_vcf import fastq_to_vcf, make_fastq
from snapperdb.snpdb import vcf_to_db, make_snpdb, get_the_snps, update_distance_matrix,\
                            update_clusters, add_ref_cluster


def setup_logging(args):
  
  """
  Function to set up logging


  Parameters
  ----------
  args : argpase object
      object containing userpassed arguments

  Returns
  -------
  logger : 
      logger object


  """

        
  logger = logging.getLogger('snapperdb')
  args.now = str(datetime.datetime.now()).split('.')[0].replace(' ', '_').replace(':', '.')
  logging.basicConfig(filename='%s/%s.snapperdb.log' % (args.log_dir, args.now), level=logging.DEBUG,
                      format='%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s')
  return logger

# -------------------------------------------------------------------------------------------------


def run_command(args):
    config_dict = parse_config(args)
    setup_logging(args)

    if args.command == 'fastq_to_db':
        logger = logging.getLogger('snapperdb.fastq_to_db')
        args.fastqs = sorted(args.fastqs)
        logger.info('PARAMS: config = %s; fastqs = %s' % (args.config_file, str(args.fastqs)))
        vcf = fastq_to_vcf(args, config_dict)
        vcf_to_db(args, config_dict, vcf)
    elif args.command == 'fastq_to_vcf':
        logger = logging.getLogger('snapperdb.fastq_to_vcf')
        args.fastqs = sorted(args.fastqs)
        logger.info('PARAMS: config = %s; fastqs = %s' % (args.config_file, str(args.fastqs)))
        vcf = fastq_to_vcf(args, config_dict)
    elif args.command == 'vcf_to_db':
        logger = logging.getLogger('snapperdb.vcf_to_db')
        logger.info('PARAMS: config = %s; vcf = %s' % (args.config_file, args.vcf))
        # third argument is for an instance of a vcf class, which doesnt exist in this case
        vcf_to_db(args, config_dict, None)

    elif args.command == 'update_distance_matrix':
        logger = logging.getLogger('snapperdb.update_distance_matrix')
        logger.info('PARAMS: config = %s;' % (args.config_file))
        update_distance_matrix(config_dict, args)

    elif args.command == 'update_clusters':
        logger = logging.getLogger('snapperdb.update_clusters')
        logger.info('PARAMS: config = %s' % args.config_file)
        update_clusters(config_dict)

    elif args.command == 'make_snpdb':
        logger = logging.getLogger('snapperdb.make_snpdb')
        logger.info('PARAMS: config = %s' % args.config_file)
        make_snpdb(config_dict)
        #map reference against itself wgsim some FASTQs if not provided
        vcf = make_fastq(args, config_dict)
        vcf_to_db(args, config_dict, vcf)
        #add reference genome to strain_clusters
        add_ref_cluster(args,config_dict)

    elif args.command == 'get_the_snps':
        logger = logging.getLogger('snapperdb.get_the_snps')
        logger.info(
            'PARAMS: config = %s; list = %s; snp_cutoff = %s' % (args.config_file, args.strain_list, args.snp_co))
        if args.out == None:
            args.out = '%s.%s.%s' % (
            str(datetime.datetime.now()).split(' ')[0], config_dict['snpdb_name'], args.strain_list)
        if args.rec_file != 'N' and args.gubbins_rec_file != None:
            logger.error('Please only pass one of rec_file or gubbins_rec_file.')
            print 'Please only pass one of rec_file or gubbins_rec_file.'
            sys.exit()
        get_the_snps(args, config_dict)

# -------------------------------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(prog='snapperdb')
    parser.add_argument("-v", "--version", help="Installed snapperdb version", action="version",
                        version="%(prog)s " + str(__version__))

    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

    parser_fastq_to_db = subparsers.add_parser('fastq_to_db', help='Takes fastqs and a config file, produces vcf and '
                                                                   'then inserts variants into the snpdb specified '
                                                                   'in the config')
    parser_fastq_to_db.add_argument('fastqs', metavar='FASTQ file(s)', nargs='+', help='At least one fastq file')
    parser_fastq_to_db.add_argument('-c', dest='config_file', metavar='Config file', required=True,
                                    help='The name of a config file in the user_configs directory (not the full path)')
    parser_fastq_to_db.add_argument('-g', dest='log_dir', default=os.path.join(os.path.expanduser('~'), 'logs'),
                                     help='Where do you want the logs written to? Will default to /path/to/fastq/logs')    
    parser_fastq_to_vcf = subparsers.add_parser('fastq_to_vcf', help='Takes fastqs and a config file and produces a '
                                                                     'vcf and '
                                                                     'serialised SNPs and ignored positions')
    parser_fastq_to_vcf.add_argument('fastqs', metavar='FASTQ file(s)', nargs='+', help='At least one fastq file')
    parser_fastq_to_vcf.add_argument('-c', dest='config_file', metavar='Config file', required=True,
                                     help='The name of a config file in the user_configs directory (not the full path)')
    parser_fastq_to_vcf.add_argument('-g', dest='log_dir', default=os.path.join(os.path.expanduser('~'), 'logs'),
                                     help='Where do you want the logs written to? Will default to /path/to/fastq/logs')
    parser_vcf_to_db = subparsers.add_parser('vcf_to_db', help='Takes a vcf and a config file, parses the vcf and then '
                                                               'adds to '
                                                               'snpdb specified in the config file')
    parser_vcf_to_db.add_argument('vcf', metavar='VCF file', nargs='+', help='A vcf file (generated using '
                                                                             'emit_all_positions?)')
    parser_vcf_to_db.add_argument('-c', dest='config_file', metavar='Config file', required=True,
                                  help='The name of a config file in the user_configs directory (or the contig path)')
    parser_vcf_to_db.add_argument('-g', dest='log_dir', default=os.path.join(os.path.expanduser('~'), 'logs'),
                                  help='Where do you want the logs written to? Will default to /path/to/fastq/logs')
    parser_make_snpdb = subparsers.add_parser('make_snpdb', help='Takes a config and makes a snpdb')
    parser_make_snpdb.add_argument('-c', dest='config_file', metavar='Config file', required=True,
                                   help='The name of a config file in the user_configs directory (not the full path)')
    parser_make_snpdb.add_argument('-g', dest = 'log_dir', default = os.path.join(os.path.expanduser('~'), 'logs'),
                                     help='Where do you want the logs written to? Will default to /path/to/fastq/logs')
    parser_make_snpdb.add_argument('-f', dest = 'fastqs', default = [], nargs='+',
                                     help='This should be ignored - I DONT KNOW HOW TO ADD FASTQ TO NAMESPACE')
    parser_update_distance_matrix = subparsers.add_parser('update_distance_matrix',
                                                          help='Takes a config and updates the distance matrix in the '
                                                               'specified snpdb')
    parser_update_distance_matrix.add_argument('-c', dest='config_file', metavar='Config file', required=True,
                                               help='The name of a config file in the user_configs directory '
                                                    '(not the full path)')
    parser_update_distance_matrix.add_argument('-g', dest='log_dir',
                                               default=os.path.join(os.path.expanduser('~'), 'logs'),
                                               help='Where do you want the logs written to? Will default to a '
                                                    '/user/home/logs')
    parser_update_distance_matrix.add_argument('-n', dest='now', default='N',
                                               help='')

    parser_get_the_snps = subparsers.add_parser('get_the_snps',
                                                help='Takes a config file, a list, and a bunch of other flags and '
                                                     'provides you with snps and more')
    parser_get_the_snps.add_argument('-c', dest='config_file', metavar='Config file', required=True,
                                     help='The name of a config file in the user_configs directory (not the full path)')
    parser_get_the_snps.add_argument('-l', dest='strain_list', required=True)

    parser_get_the_snps.add_argument('-m', dest='snp_co', required=True,
                                     help='SNP cut off (has to be integer), strains more than this number of '
                                          'SNPs from the reference will be excluded from the analysis. '
                                          'A sensible starting point is 3000')
    parser_get_the_snps.add_argument('-o', dest='out', help='Prefix for output, will default to '
                                                            '<date>.<snpdb_name>_vs_<list_name>')
    parser_get_the_snps.add_argument('-a', dest='alignment_type',
                                     help='Alignment type (W=whole Consensus, A=Soft Core, C=Core), default is alignments where at least 80%% of the samples are A/C/T/G at each position (A:80)',
                                     default='A:80')
    parser_get_the_snps.add_argument('-x', dest='mat_flag',
                                     help='Would you like a pairwise distance matrix? Y/N (default = N)', default='N')
    parser_get_the_snps.add_argument('-v', dest='var_flag', help='Would you like a more detailed list of the variant '
                                                                 'attributes? Y/N (default = N)')
    parser_get_the_snps.add_argument('-r', dest='ref_flag',
                                     help='Would you like the reference genome in the alignment? Y/N ('
                                          'default = N)', default='N')
    parser_get_the_snps.add_argument('-n', dest='rec_file',
                                     help='File with list of positions to ignore due to expected '
                                          'recombination', default='N')
    parser_get_the_snps.add_argument('-ng', dest = 'gubbins_rec_file', help = 'The gff file produced by gubbins with '
                                    'recombinant positions in it.')
    parser_get_the_snps.add_argument('-b', dest='back_flag',
                                     help='Would you like a background cluster level? SNP cluster level '
                                          'from which to take one representative', default='N')
    parser_get_the_snps.add_argument('-g', dest='log_dir', default=os.path.join(os.path.expanduser('~'), 'logs'),
                                     help='Where do you want the logs written to? Will default to a /user/home/logs')
    parser_update_clusters = subparsers.add_parser('update_clusters',
                                                   help='Given a config file, updates the SNP clustering '
                                                        'associated with the SNPdb specified in the config.')
    parser_update_clusters.add_argument('-c', dest='config_file', metavar='Config file', required=True,
                                        help='The name of a config '
                                             'file in the user_configs directory (not the full path)')
    parser_update_clusters.add_argument('-g', dest='log_dir', default=os.path.join(os.path.expanduser('~'), 'logs'),
                                        help='Where do you want the logs written to? Will default to a /user/home/logs')


    args = parser.parse_args()
    print args
    run_command(args)

    print "### Completed "+str(datetime.datetime.now())

if __name__ == "__main__":
    main()
