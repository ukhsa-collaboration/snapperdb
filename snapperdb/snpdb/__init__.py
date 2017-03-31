__author__ = 'gidis'

from datetime import datetime
import inspect
import os
import re
import sys
import logging
import psycopg2, psycopg2.extras
from snpdb import SNPdb
import snapperdb
from snapperdb.gbru_vcf import Vcf
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
    #annotate vars
    logger.info('Annotating new variants')

    snpdb.snpdb_annotate_vars(vcf)

# -------------------------------------------------------------------------------------------------

def add_ref_cluster(args, config_dict):
    #set up loggging
    logger = logging.getLogger('snapperdb.snpdb.add_ref_cluster')
    logger.info('Initialising SNPdb class')

    #create snpdb class
    snpdb = SNPdb(config_dict)

    #parse config into snpdb object
    logger.info('Parsing config dict')
    snpdb.parse_config_dict(config_dict)
    
    #connect to snpdb postgres
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)

    snpdb.add_cluster()

# -------------------------------------------------------------------------------------------------

def make_snpdb(config_dict):
    snpdb = SNPdb(config_dict)
    snpdb._connect_to_snpdb()
    snpdb.make_snpdb()

# -------------------------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------------------------

def read_rec_file_mc(rec_file):
    #read recombination file - tab delineated with contig
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

# -------------------------------------------------------------------------------------------------

def create_contig_index_for_consensus_genome(reference_genome):
    ## first, need to concatenate the reference genome in the order of 
    ## alphanumerically sorted contig names
    concatenated_ref_genome = ''
    contig_names = sorted(reference_genome.keys())
    for c in contig_names:
        concatenated_ref_genome += ''.join(reference_genome[c])
    ## contig dict is going to be {(contig start in concat ref genome, contig stop in concat ref genome):contig_name}
    contig_index = {}
    ## start at 0
    i = 0
    ## need to iterate through the contig names in the sorted order, because this is the order that get_the_snps will go through them in
    for contig in contig_names:
        ## use i and j to move through ref genome
        ## we minus 1 here because of the python thing of not counting the last 'fence post' in the index, but counting it in len(). I think. Seems to work.
        j = i + len(reference_genome[contig]) - 1
        contig_index[(i, j)] = contig
        i += len(reference_genome[contig])
        
    pprint.pprint(contig_index)
    return contig_index

# -------------------------------------------------------------------------------------------------

def make_recomb_dict_from_gubbins(recombinant_sections, contig_index):
    rec_dict = {}
    for rs in recombinant_sections:
        for contig in contig_index:
            if contig[0] <= rs[0] <= contig[1]:
                if contig[0] <= rs[1] <= contig[1]:
                    recomb_start = (rs[0] - contig[0])
                    recomb_stop = (rs[1] - contig[0])
                    if contig_index[contig] in rec_dict:
                        rec_dict[contig_index[contig]] += range(recomb_start, recomb_stop)
                    else:
                        rec_dict[contig_index[contig]] = []
                        rec_dict[contig_index[contig]] += range(recomb_start, recomb_stop)
                else:
                    print 'The recombination section is not on one contig, this may not be a real recombination event, this is not currently being exlcuded.', rs, contig, contig_index[contig]
    for contig in rec_dict:
        rec_dict[contig] = set(rec_dict[contig])
    return rec_dict

# -------------------------------------------------------------------------------------------------

def read_rec_file_mc_gubbins(gubbins_rec_file, reference_genome):
    contig_index = create_contig_index_for_consensus_genome(reference_genome)
    try:
        openfile = open(gubbins_rec_file, 'r')
    except IOError:
        print gubbins_rec_file + ' not found ...'
        sys.exit()
    recombinant_sections = []
    for line in openfile.readlines():
        if not line.startswith('##'):
            split_line = line.strip().split('\t')
            recombinant_sections.append((int(split_line[3]), int(split_line[4])))
            
    rec_dict = make_recomb_dict_from_gubbins(recombinant_sections, contig_index)
    return rec_dict
    
 # -------------------------------------------------------------------------------------------------
   

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
    

    #query snadb    
    snpdb.parse_args_for_get_the_snps_mc(args, strain_list, ref_seq, snpdb.reference_genome)
    
    snpdb.print_fasta_mc(args, rec_dict)

    #print matrix
    if args.mat_flag == 'Y':
        snpdb.print_matrix(args.out)
    # print variant list    
    if args.var_flag == 'Y':
        logger.info('Printing variants')
        snpdb.print_vars_mc(args,rec_dict)

# -------------------------------------------------------------------------------------------------

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
        print '### Launching serial update_distance_matrix ' + str(datetime.now())
        snpdb.check_matrix(strain_list, update_strain)
        snpdb.update_clusters()
    else:
        print '### Nothing to update ' + str(datetime.now())
# -------------------------------------------------------------------------------------------------


def update_clusters(config_dict):
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    snpdb.update_clusters()
