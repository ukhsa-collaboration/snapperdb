__author__ = 'gidis'

from datetime import datetime
import inspect
import os
import pickle
import re
import sys
import logging
import psycopg2, psycopg2.extras
from snpdb import SNPdb
import snapperdb
from snapperdb.gbru_vcf import Vcf
import glob
import pprint


def vcf_to_db(args, config_dict, vcf):
    #set up loggging
    logger = logging.getLogger('snapperdb.snpdb.vcf_to_db')
    logger.info('Initialising SNPdb class')

    #create snpdb class
    snpdb = SNPdb(config_dict)

    #parse config into snpdb object
    logger.info('Parsing config dict')
    snpdb.parse_config_dict(config_dict)
    
    #connect to snpdb postgres
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)

    #check stack?
    if inspect.stack()[0][3] == 'fastq_to_db':
        # fastq_to_db we will alread have a vcf object to work wih
        logger.info('You are running fastq_to_db.')
        
    elif inspect.stack()[0][3] == 'vcf_to_db':
        ## there is no existing vcf class here, but there will definitely be a vcf
        logger.info('You are running vcf_to_db. Initialising Vcf class.')
        vcf = Vcf()
        logger.info('Making SNPdb variables and output files')
        #set up variables
        snpdb.define_class_variables_and_make_output_files(args, vcf)
    
    #read vcf
    vcf.read_multi_contig_vcf()
    logger.info('Uploading to SNPdb')
    #upload vcf
    snpdb.snpdb_upload(vcf)


def make_snpdb(config_dict):
    snpdb = SNPdb(config_dict)
    snpdb._connect_to_snpdb()
    snpdb.make_snpdb()

def read_file(file_name):
    #read list of strains for get the snps
    try:
        openfile = open(file_name, 'r')
    except:
        print file_name + " not found ... "
        sys.exit()
    strain_list = []
    for line in openfile:
        strain_list.append(line.strip())
    return strain_list


def read_multi_contig_fasta(ref):
    try:
        openfile = open(ref, 'r')
    except:
        print (ref + " not found ... ")
        sys.exit()
    ref_seq = {}
    contig = ""
    for line in openfile:
        matchObj = re.search('>', line)
        if matchObj is None:
            for n in line.strip():
                ref_seq[contig[0]].append(n)
        else:
            contig = line[1:].strip().split()
            ref_seq[contig[0]] = []
    return ref_seq


def read_rec_file_mc(rec_file):
    #read recombination file - tab delineated with reference genome
    try:
        openfile = open(rec_file, 'r')
    except:
        print (rec_file + " not found ... ")
        sys.exit()
    
    rec_dict = {}
    for line in openfile:
        split_line = line.strip().split('\t')
        rec_range = range(int(split_line[1]) - 1, (int(split_line[2]) - 1))
        if split_line[0] in rec_dict:
            rec_dict[split_line[0]] = set(rec_dict[split_line[0]]) | set(rec_range)
        else:
            rec_dict[split_line[0]] = set(rec_range)
    
    return rec_dict

def read_rec_file_mc_gubbins(gubbins_rec_file, reference_genome):
    ## first, need to concatenate the reference genome in the order of 
    ## alphanumerically sorted contig names
    concatenated_ref_genome = ''
    contig_names = sorted(reference_genome.keys())
    for c in contig_names:
        # print len(reference_genome[c])
        concatenated_ref_genome += ''.join(reference_genome[c])
    print len(concatenated_ref_genome)

    ## contig dict is going to be {(contig start in concat ref genome, contig stop in concat ref genome):contig_name}
    contig_dict = {}
    ## start at 0
    i = 0
    for contig in reference_genome:
        ## use i and j to move through ref genome
        ## we minus 1 here because of the python thing of not counting the last 'fence post' in the index, but counting it in len()
        j = i + len(reference_genome[contig]) - 1
        contig_dict[(i, j)] = contig
        i += len(reference_genome[contig])
        
    pprint.pprint(contig_dict)
    try:
        openfile = open(gubbins_rec_file, 'r')
    except IOError:
        print gubbins_rec_file + ' not found ...'
        sys.exit()
    

    

def get_the_snps(args, config_dict):
    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.get_the_snps')
    logger.info('Inititialising SnpDB Class')
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)
    #read strainlist
    strain_list = read_file(args.strain_list)

    #connect to postgresdb
    snpdb._connect_to_snpdb()
    #get reference genome path
    ref_seq_file = os.path.join(snapperdb.__ref_genome_dir__, snpdb.reference_genome + '.fa')
    #read the reference fasta
    ref_seq = read_multi_contig_fasta(ref_seq_file)
    #add reference genome to strain_list
    strain_list.append(snpdb.reference_genome)

    #if recombination flag set
    if args.rec_file != 'N':
        logger.info('Reading recombination list')
        rec_dict = read_rec_file_mc(args.rec_file)
    elif args.gubbins_rec_file != None:
        logger.info('Reading gubbins recombination list')
        rec_dict = read_rec_file_mc_gubbins(args.gubbins_rec_file, ref_seq)
    else:
        #should we set this as none
        rec_dict = {}
    
    # sys.exit()

    #query snadb    
    snpdb.parse_args_for_get_the_snps_mc(args, strain_list, ref_seq, snpdb.reference_genome)
    #print fasta
    snpdb.print_fasta_mc(args, rec_dict)

    #print matrix
    if args.mat_flag == 'Y':
        snpdb.print_matrix(args.out)
    # print variant list    
    if args.var_flag == 'Y':
        logger.info('Printing variants')
        snpdb.print_vars_mc(args,rec_dict)




def update_distance_matrix(config_dict, args):
    logger = logging.getLogger('snapperdb.snpdb.update_distance_matrix')
    logger.info('Inititialising SnpDB Class')
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    logger.info('Getting strains')
    strain_list, update_strain = snpdb.get_strains()
    # # get_all_good_ids from snpdb2 takes a snp cutoff as well, here, we don't have a SNP cutoff so we set it arbitrarily high.
    snp_co = '1000000'
    if update_strain:
        print "###  Populating distance matrix: " + str(datetime.now())
        snpdb.parse_args_for_update_matrix(snp_co, strain_list)
        if args.hpc == 'N':
            print '### Launching serial update_distance_matrix ' + str(datetime.now())
            snpdb.check_matrix(strain_list, update_strain)
            snpdb.update_clusters()
        else:
            try:
                print '### Launching parallel update_distance_matrix ' + str(datetime.now())
                args.hpc = int(args.hpc)
                short_strain_list = set(strain_list) - set(update_strain)
                snpdb.write_qsubs_to_check_matrix(args, strain_list, short_strain_list, update_strain, config_dict['snpdb_name'])
                # # on cluster version this will have to be subject to a qsub hold - no it wont, can just run on headnode
                snpdb.check_matrix(update_strain, update_strain)
            except ValueError as e:
                print '\n#### Error ####'
                print e, '-m has to be an integer'
    else:
        print '### Nothing to update ' + str(datetime.now())

def qsub_to_check_matrix(config_dict, args):
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    snp_co = '1000000'
    strain_list = []
    with open(args.strain_list) as fi:
        for x in fi.readlines():
            strain_list.append(x.strip())
    short_strain_list = []
    with open(args.short_strain_list) as fi:
        for x in fi.readlines():
            short_strain_list.append(x.strip())
    update_strain = []
    with open(args.update_list) as fi:
        for x in fi.readlines():
            update_strain.append(x.strip())
    snpdb.parse_args_for_update_matrix(snp_co, strain_list)
    snpdb.check_matrix(short_strain_list, update_strain)

    # # need to clean up as otherwise the glob
    os.system('rm -f {0}'.format(args.strain_list))
    direc, name = os.path.split(args.strain_list)
    list_number = name.split('_')[-1]
    shell_script = '{0}/update_matrix_{1}.sh'.format(direc, list_number)
    os.system('rm -f {0}'.format(shell_script))

def update_clusters(config_dict):
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    snpdb.update_clusters()

def get_variants_of_interest(config_dict, args):
    background_list = read_file(args.background_list)
    of_interest_list = read_file(args.of_interest_list)


    '''
    To do
    1. for each list, get the good quality variants for each isolate into {strain:[good, vars], ...}
    2.

    '''

def upload_indels(config_dict, args):
    snpdb = SNPdb(config_dict)
    snpdb._connect_to_snpdb()
    vcf = Vcf()
    vcf.parse_config_dict(config_dict)
    snpdb.define_class_variables_and_make_output_files_indels(args, vcf)
    snpdb.add_indels_to_snpdb(vcf)



    '''
    1. parse vcf for indels - have different function based on parse_vcf for now
        a. check that length of ref and alt are different
        b. then apply normal quality filters - dp, ad, gq

    2. get all indels from indels table in snpdb
        a. run equivalent of snpdb.add_to_snpdb for indels

    science - need to check how often indels start at the same or similar (within a few positions) positions

    '''
    pass


