__author__ = 'flashton'

import os
import errno


# def mkdir_p(path):
#     try:
#         os.makedirs(path)
#     except OSError as exc:
#         if exc.errno == errno.EEXIST and os.path.isdir(path):
#             pass
#         else:
#             raise
#

# def main():

    # ## here, we call the vcf class and assign to variable bam as we are going to be doing more bam-y things
    # fastq_bam_vcf = Vcf()
    # fastq_bam_vcf.path_to_config = path_to_config
    # fastq_bam_vcf.sample_name = os.path.basename(opts.fastq_1).split(os.extsep)[0]
    # fastq_bam_vcf.make_sorted_bam(opts)
    # fastq_bam_vcf.make_vcf(opts)


if __name__ == '__main__':
    '''
    This will be deprecated, as should only access through SnapperDB_main.py
    '''

    # parser = argparse.ArgumentParser(description='The function of this script is to take a config file and path to vcf and '
    #                                              'upload the SNPs and ignored positions to SNPdb',
    #                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('--config', action='store', help='config file', required=True)
    # parser.add_argument('--fastq_1', action='store', help='fastq_1', required=True)
    # parser.add_argument('--fastq_2', action='store', help='fastq_2', required=True)
    # # not expected from user, will be derived from the fastq path. just added here because putting it in opts makes it easy to
    # #  pass around
    # parser.add_argument('--out_dir', action='store', help='out_dir (e.g. tmp)')

    # opts = parser.parse_args()
    # opts.config = os.path.basename(opts.config)
    # path_to_config = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'user_configs', opts.config))
    # opts.out_dir = os.path.join(os.path.dirname(opts.fastq_1), 'snpdb', 'tmp')
    # mkdir_p(opts.out_dir)

    # fastq_to_vcf(opts, path_to_config)