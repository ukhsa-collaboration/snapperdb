__author__ = 'flashton'

import sys
import os
import argparse

from snpdb_modules import SNPdb


"""
User to run this first, takes config file

"""


useage = 'Fill this in with something useful\n'


def main():

    parser = argparse.ArgumentParser(description='The function of this script is to take a config file and path to vcf and '
                                                 'upload the SNPs and ignored positions to SNPdb',
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--config', action = 'store', help = 'config file')

    ## optionally, can pass all the stuff from the config file separately
    parser.add_argument('--snpdb_name', action = 'store', help = 'snpdb_name')
    parser.add_argument('--reference_genome', action = 'store', help = 'reference_genome')
    parser.add_argument('--pg_uname', action = 'store', help = 'pg_uname')
    parser.add_argument('--pg_pword', action = 'store', help = 'pg_pword')
    parser.add_argument('--pg_host', action = 'store', help = 'pg_host')
    opts = parser.parse_args()


    ## this ensures that if the user passes the absolute or relative path, that path_to_config works nicely



    if opts.config:
        opts.config = os.path.basename(opts.config)
        path_to_config = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'snpdb_configs', opts.config)
        db = SNPdb(config=path_to_config, make_db=True)
    elif opts.snpdb_name:
        db = SNPdb(snpdb_name=opts.snpdb_name, reference_genome=opts.reference_genome, pg_uname=opts.pg_uname,
                   pg_pword=opts.pg_pword, pg_host=opts.pg_host, make_db=True)
    else:
        sys.stderr.write(useage)


if __name__ == '__main__':
    main()