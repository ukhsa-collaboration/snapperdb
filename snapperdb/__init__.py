__author__ = 'flashton'

import os
import sys
import subprocess

## in order to include the commit number in the version, need to
## store the current working direcotory
#original_wd = os.getcwd()
## change wd to the parent script directory e.g. snapperdb
#os.chdir(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
## the version is digits + the output of the git rev-parse function
__version__ = "v0.1.1 " # + subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
## change wd back to original_wd
#os.chdir(original_wd)


__config_dir__ = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'user_configs')
__ref_genome_dir__ = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'reference_genomes')


def parse_config(args):
    path_to_config = os.path.join(__config_dir__, args.config_file) # - deprecating this as need to handle not having configs
    ## in the user_configs that can be found relative to the script dir.
    ## however, can still use just the name of the file as user_configs in the module dir is in the path when you load the module
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
                    config_dict['depth_cutoff'] = float(line.strip().split()[-1])
                if line.startswith('mq_cutoff'):
                    config_dict['mq_cutoff'] = float(line.strip().split()[-1])
                if line.startswith('ad_cutoff'):
                    config_dict['ad_cutoff'] = float(line.strip().split(' ')[-1])
                if line.startswith('average_depth_cutoff'):
                    config_dict['average_depth_cutoff'] = float(line.strip().split()[-1])
    except IOError:
        print 'Cannot find {0}'.format(path_to_config)
        sys.exit()

    return config_dict

    # # don't need to check the below here necassarily, probably checked later.
    # ref_genome_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'reference_genomes')
    # if os.path.exists(ref_genome_dir):
    #     pass
    # else:
    #     sys.stderr.write('Ref genome dir not found\n')
    #     sys.exit()
