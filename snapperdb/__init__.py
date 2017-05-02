__author__ = 'gidis'

import os
import sys


__version__ = "v1.0 "  


__config_dir__ = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'user_configs')

# in the brave new world we no longer keep stuff in relative paths to the code
# the location of the reference genomes is set in the module for phe/gastro_fastq_to_vcf_config
# keep the except branch for compatibility - ulf 1Nov2016
try:
    __ref_genome_dir__ = os.path.join(os.environ['GASTROSNAPPER_REFPATH'], "1.0")
except KeyError:
    __ref_genome_dir__ = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'reference_genomes')


def parse_config(args):
    if os.path.isabs(args.config_file) == False:
        path_to_config = os.path.join(__config_dir__, args.config_file) # - deprecating this as need to handle not having configs
    else:
        path_to_config = args.config_file
    # in the user_configs that can be found relative to the script dir.
    # however, can still use just the name of the file as user_configs in the module dir is in the
    # path when you load the module

    #add mapper, variant top config - TD 270317

    config_dict = {}
    try:
        with open(path_to_config, 'r') as fi:
            for line in fi.readlines():
                if line.startswith('reference_genome'):
                    config_dict['reference_genome'] = line.strip().split()[-1]
                if line.startswith('snpdb_name'):
                    config_dict['snpdb_name'] = line.strip().split()[-1]
                if line.startswith('pg_uname'):
                    config_dict['pg_uname'] = line.strip().split()[-1]
                if line.startswith('pg_pword'):
                    config_dict['pg_pword'] = line.strip().split()[-1]
                if line.startswith('pg_host'):
                    config_dict['pg_host'] = line.strip().split()[-1]
                if line.startswith('depth_cutoff'):
                    config_dict['depth_cutoff'] = int(line.strip().split()[-1])
                if line.startswith('mq_cutoff'):
                    config_dict['mq_cutoff'] = int(line.strip().split()[-1])
                if line.startswith('ad_cutoff'):
                    config_dict['ad_cutoff'] = float(line.strip().split(' ')[-1])
                if line.startswith('average_depth_cutoff'):
                    config_dict['average_depth_cutoff'] = int(line.strip().split()[-1])
                if line.startswith('mapper'):
                    config_dict['mapper'] = line.strip().split()[-1]
                if line.startswith('variant_caller'):
                    config_dict['variant_caller'] = line.strip().split()[-1]           

    except IOError:
        print 'Cannot find {0}'.format(path_to_config)
        sys.exit()
    return config_dict

# -------------------------------------------------------------------------------------------------


def parse_ancillary_info(infile):
    res_dict = {}
    with open(infile, 'r') as fi:
        for line in fi.readlines():
            split = line.strip().split('\t')
            res_dict[split[0]] = split[1]
    return res_dict
