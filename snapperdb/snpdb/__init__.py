__author__ = 'gidis'

import datetime
import inspect
import os
import re
import subprocess
import sys
import json
import logging
import psycopg2, psycopg2.extras
import snapperdb
from snpdb import SNPdb
from snapperdb.gbru_vcf import Vcf
from snapperdb import parse_config
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
    snpdb.snpdb_upload(vcf,args)
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
        print "### Strain list '" + file_name + "' not found ... "
        print "### Exiting "+str(datetime.datetime.now())
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
        print "### Reference genome "+ ref + " not found ... "
        print "### Exiting "+str(datetime.datetime.now())
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

def read_rec_file_mc(rec_file, rec_dict):
    #read recombination file - tab delineated with contig
    try:
        openfile = open(rec_file, 'r')
    except:
        print "### Recombination list "+ rec_file + " not found ... "
        print "### Exiting "+str(datetime.datetime.now())
        sys.exit()

    for line in openfile:
        split_line = line.strip().split('\t')
        if split_line[0] in rec_dict:
            try:
                rec_dict[split_line[0]] += range(int(split_line[1]), (int(split_line[2])))
            except:
                print 'Error parsing ignored position list'
                sys.exit()
        else:
            rec_dict[split_line[0]] = []
            try:
                rec_dict[split_line[0]] += range(int(split_line[1]), (int(split_line[2])))
            except:
                print 'Error parsing ignored position list'
                sys.exit()

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

def make_recomb_dict_from_gubbins(recombinant_sections, contig_index, rec_dict):
    for rs in recombinant_sections:
        for contig in contig_index:
            if contig[0] <= rs[0] <= contig[1]:
                if contig[0] <= rs[1] <= contig[1]:
                    recomb_start = (rs[0] - contig[0])
                    recomb_stop = (rs[1] - contig[0])
                    if contig_index[contig] in rec_dict:
                        try:
                            rec_dict[contig_index[contig]] += range(recomb_start, recomb_stop)
                        except:
                            print 'Error parsing GFF'
                            sys.exit()

                    else:
                        rec_dict[contig_index[contig]] = []
                        try:
                            rec_dict[contig_index[contig]] += range(recomb_start, recomb_stop)
                        except:
                            print 'Error parsing GFF'
                            sys.exit()
                  
                else:
                    print 'The recombination section is not on one contig, this may not be a real recombination event, this is not currently being exlcuded.', rs, contig, contig_index[contig]
    for contig in rec_dict:
        rec_dict[contig] = set(rec_dict[contig])
    return rec_dict

# -------------------------------------------------------------------------------------------------

def read_rec_file_mc_gubbins(gubbins_rec_file, reference_genome, rec_dict):
    contig_index = create_contig_index_for_consensus_genome(reference_genome)
    try:
        openfile = open(gubbins_rec_file, 'r')
    except IOError:
        print "### Gubbins GFF file "+ gubbins_rec_file + " not found ... "
        print "### Exiting "+str(datetime.datetime.now())
        sys.exit()
    recombinant_sections = []
    for line in openfile.readlines():
        if not line.startswith('##'):
            split_line = line.strip().split('\t')
            recombinant_sections.append((int(split_line[3]), int(split_line[4])))

    rec_dict = make_recomb_dict_from_gubbins(recombinant_sections, contig_index, rec_dict)
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
    rec_dict = {}

    if args.rec_file != 'N':
        logger.info('Reading recombination list')
        rec_dict = read_rec_file_mc(args.rec_file, rec_dict)
    if args.gubbins_rec_file != None:
        logger.info('Reading gubbins recombination list')
        rec_dict = read_rec_file_mc_gubbins(args.gubbins_rec_file, ref_seq, rec_dict)


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


def export_json(args, config_dict):
    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.export_json')
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


    snpdb.parse_args_for_export(args, strain_list, ref_seq)


  # -------------------------------------------------------------------------------------------------

def ignore_isolate(args,config_dict):

    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.export_json')
    logger.info('Inititialising SnpDB Class')
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)

    #connect to postgresdb
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)  

    #remove isolate
    snpdb.remove_isolate(args.ig_strain)



  # -------------------------------------------------------------------------------------------------

def get_strains(args,config_dict):

    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.get_strains')
    logger.info('Inititialising SnpDB Class')
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)

    #connect to postgresdb
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)  

    #remove isolate
    snpdb.get_strain_list(args.thresh)


  # -------------------------------------------------------------------------------------------------


def accept_outlier(args,config_dict):

    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.export_json')
    logger.info('Inititialising SnpDB Class')
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)

    #connect to postgresdb
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)  

    #remove isolate
    snpdb.zscore_exception(args.out_strain)


  # -------------------------------------------------------------------------------------------------

def import_json(args):
    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.import_json')
    json_path= untar_file(args.json_file)
    json_dict = {}

    #import json
    try:
        with open (json_path) as json_data:
            json_dict = json.load(json_data)
        json_data.close()
    except IOError:
        print "Issue with JSON file. Exiting."
        exit()

    #parse config
    args.config_file = json_dict['config_file']
    config_dict = parse_config(args)
    
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)
    #get reference genome path
    ref_seq_file = os.path.join(snapperdb.__ref_genome_dir__, snpdb.reference_genome + '.fa')
    #read the reference fasta
    ref_seq = read_multi_contig_fasta(ref_seq_file)

    #create VCF class
    vcf = Vcf()
    vcf.parse_json_dict(json_dict, ref_seq)
    vcf.depth_average = json_dict['strain_stats']
    vcf.sample_name = json_dict['sample']

    logger.info('Uploading to SNPdb')
    #upload vcf
    #connect to snpdb postgres
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)    
    
    if args.write_flag == 'W':

        snpdb.snpdb_upload(vcf,args)
        #annotate vars
        logger.info('Annotating new variants')

        snpdb.snpdb_annotate_vars(vcf)
    elif args.write_flag == 'R':
        print "under development"
        #snpdb.snpdb_query(vcf,args)

# -------------------------------------------------------------------------------------------------

def untar_file(file_path):
    if not file_path.endswith(".tar.gz"):
        print "Not a Tar file. Exiting."
        exit()
    bash_command = "tar -xzvf " + file_path
    flag = subprocess.call(bash_command, shell = True)
    if flag != 0:
        print "Issue with Tar file. Exiting."
        exit()
    new_path = file_path[0:-7]
    return new_path

# -------------------------------------------------------------------------------------------------


def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i + n]

def update_distance_matrix(config_dict, args):
    logger = logging.getLogger('snapperdb.snpdb.update_distance_matrix')
    logger.info('Inititialising SnpDB Class')
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    logger.info('Getting strains')
    strain_list, update_strain, all_strains = snpdb.get_strains()


    # # get_all_good_ids from snpdb2 takes a snp cutoff as well, here, we don't have a SNP cutoff so we set it arbitrarily high.
    snp_co = '1000000'
    if all_strains or len(update_strain) > 1:
        if update_strain:
            print "###  Populating distance matrix: " + str(datetime.datetime.now())
            snpdb.parse_args_for_update_matrix(snp_co, strain_list)
            if args.hpc == 'N':
                print '### Launching serial update_distance_matrix ' + str(datetime.datetime.now())
                snpdb.check_matrix(strain_list, update_strain)
                snpdb.update_clusters()
            else:
                print '### Launching parallel update_distance_matrix ' +str(datetime.datetime.now())
                present_stains = list(set(strain_list) - set(update_strain))
                for idx, one_strain in enumerate(chunks(list(update_strain), int(args.hpc))):
                    snpdb.write_qsubs_to_check_matrix(args, idx, one_strain, present_stains, config_dict['snpdb_name'])
                snpdb.check_matrix(update_strain, update_strain)
        else:
            print '### Nothing to update ' + str(datetime.datetime.now())
    else:
        print '### Nothing to update ' + str(datetime.datetime.now())
# -------------------------------------------------------------------------------------------------

def qsub_to_check_matrix(config_dict, args):

    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    snp_co = '1000000'
    added_list = []
    with open(args.added_list) as fi:
        for x in fi.readlines():
            added_list.append(x.strip())
    present_strains = []
    with open(args.present_strains) as fi:
        for x in fi.readlines():
            present_strains.append(x.strip())

    strain_list = list(set(present_strains) | set(added_list))
    snpdb.parse_args_for_update_matrix(snp_co, strain_list)
    snpdb.check_matrix(strain_list,added_list)


# -------------------------------------------------------------------------------------------------


def update_clusters(config_dict):
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    snpdb.update_clusters()