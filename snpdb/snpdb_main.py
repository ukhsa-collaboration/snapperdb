__author__ = 'flashton'

import os
import sys
import argparse

from snpdb_modules import SNPdb
import gbru_vcf_modules


'''
From command line, it takes a config file (which will have to be made by the user before hand, same config file as the
gbru_vcf), path_to_vcf
'''

def main():
    parser = argparse.ArgumentParser(description='The function of this script is to take a config file and path to vcf and '
                                                 'upload the SNPs and ignored positions to SNPdb',
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--config', action = 'store', help = 'config file')
    parser.add_argument('--vcf', action = 'store', help = 'vcf file', required = True)
    ## optionally, can pass all the stuff from the config file separately
    parser.add_argument('--snpdb_name', action = 'store', help = 'snpdb_name')
    parser.add_argument('--reference_genome', action = 'store', help = 'reference_genome')
    parser.add_argument('--pg_uname', action = 'store', help = 'pg_uname')
    parser.add_argument('--pg_pword', action = 'store', help = 'pg_pword')
    parser.add_argument('--pg_host', action = 'store', help = 'pg_host')
    parser.add_argument('--depth_cutoff', action = 'store', help = 'depth_cutoff')
    parser.add_argument('--mq_cutoff', action = 'store', help = 'mq_cutoff')
    parser.add_argument('--ad_cutoff', action = 'store', help = 'ad_cutoff')
    parser.add_argument('--bad_pos_pick', action = 'store', help = 'bad_pos_pick')
    parser.add_argument('--good_var_pick', action = 'store', help = 'good_var_pick')
    parser.add_argument('--cluster', action = 'store', help = 'cluster')

    opts = parser.parse_args()

    ## thinking about removing the option for passing the separate options on the command line

    vcf = gbru_vcf_modules.Vcf()

    if opts.config:
        ## this ensures that if the user passes the absolute or relative path, that path_to_config works nicely
        opts.config = os.path.basename(opts.config)
        path_to_config = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname( __file__ )),
                                                      'user_configs', opts.config))
        #path_to_config = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'snpdb_configs', opts.config)
        if os.path.exists(path_to_config):
            pass
        else:
            sys.stderr.write('Config file doesnt exist\n')
            sys.exit()
        db = SNPdb(config=path_to_config)
        if opts.vcf.endswith('.gz'):
            os.system('gunzip {0}'.format(opts.vcf))
            opts.vcf = os.path.splitext(opts.vcf)[0]
        ## as the cutoffs are part of the SNPdb, need to pass the SNPdb instance
        ## add a depth calculating thing to the vcf class
        vcf.read_vcf(opts.vcf, db)

    elif opts.snpdb_name:
        db = SNPdb(snpdb_name=opts.snpdb_name, reference_genome=opts.reference_genome, pg_uname=opts.pg_uname,
                   pg_pword=opts.pg_pword, pg_host=opts.pg_host, depth_cutoff=opts.depth_cutoff, mq_cutoff=opts.mq_cutoff,
                   ad_cutoff=opts.ad_cutoff)

    if opts.cluster == 'True':
        ## serialise and quit as this is being done in parallel.
        ## have job which monitors the mapping, snp calling jobs that have been launched, then does batch upload.
        ## after serial upload, need to parallise the update matrix stuff.
        ## need to have an equivalent of check_len_vcf in this as well
        vcf.pickle_variants_and_ignored_pos(opts)
        sys.exit()

    else:
        db.check_len_vcf(vcf)
        db.snpdb_upload(vcf)
        ## add in the matrix update function
        db.snpdb_conn.close()

    ## as long as you have __enter__ and __exit__ methods
    ## see http://158.119.147.111/lims/lims/blob/speedup_v2/modules/uk/gov/phe/utils/db/DbUtils.py
    #with SNPdb(config="/bollocks") :
    #    db.make_db()

if __name__ == "__main__":
    main()






