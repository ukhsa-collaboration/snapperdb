__author__ = 'gidis'

import errno
import os
import pickle
import re
import subprocess
import sys
from Bio import SeqIO

import snapperdb

class ParsedVcf:
    def __init__(self):
        self.depth = {}
        self.qual = {}
        self.hap_qual = {}
        self.hap_depth = {}
        self.hap_call = {}
        self.hap_var_count = {}
        self.filter_flag = {}
        self.var = {}
        self.ref_base = {}
        self.ref = None
        self.mixed_positions = []
        self.bad_pos = []
        self.good_var = []



class Vcf:
    def __init__(self):
        self.bad_depth = None
        self.bad_qual = None
        self.bad_var = None
        self.good_var = None
        self.bad_pos = None
        self.sample_name = None
        self.depth_average = None
        self.depth_sd = None
        self.path_to_config = None
        self.reference_genome = None
        self.ref_genome_path = None
        self.snpdb_name = None
        self.tmp_dir = None
        self.vcf_filehandle = None
        self.depth_cutoff = None
        self.mq_cutoff = None
        self.ad_cutoff = None
        self.number_mixed_positions = None
        self.mixed_positions = []
        self.rec_list = []
        self.vcf_max_pos = None
        self.ref = None
        self.contig = None
        self.parsed_vcf_container = []

    def parse_config_dict(self, config_dict):
        # # we loop through thusly in case not all these things are in the config
        for attr in config_dict:
            if attr == 'reference_genome':
                self.reference_genome = config_dict[attr]
            if attr == 'depth_cutoff':
                self.depth_cutoff = config_dict[attr]
            if attr == 'mq_cutoff':
                self.mq_cutoff = config_dict[attr]
            if attr == 'ad_cutoff':
                self.ad_cutoff = config_dict[attr]

    def mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def make_tmp_dir(self, args):
        try:
            self.tmp_dir = os.path.join(os.path.dirname(args.fastqs[0]), 'snpdb', 'tmp')
            self.mkdir_p(self.tmp_dir)
        except AttributeError:
            # # this is because we also call this from snpdb, where args doesn't contain fastqs, but a pickle.
            self.tmp_dir = os.path.join(os.path.dirname(args.vcf[0]), 'tmp')
            self.mkdir_p(self.tmp_dir)


    def read_multi_contig_vcf(self):

        try:
            os.path.exists(self.vcf_filehandle)
        except IOError:
            print self.vcf_filehandle + " not found ... "

        if self.vcf_filehandle.endswith('.gz'):
            os.system('gunzip {0}'.format(self.vcf_filehandle))
            self.vcf_filehandle = os.path.splitext(self.vcf_filehandle)[0]
        try:
            openfile = open(self.vcf_filehandle, 'r')
        except:
            print self.vcf_filehandle + " not found ... "
            sys.exit()

        oref = ''

        parsed_vcf = ""

        for line in openfile:
            if not re.match(r'^#', line):
                split_line = line.split()
                # get reference
                ref = split_line[0]

                # if the reference is new
                if ref != oref:
                    # if the reference is new and this is not the first reference and the parsed vcf object
                    if oref != '':
                        self.parsed_vcf_container.append(parsed_vcf)
                    #create a parsed vcf object    
                    parsed_vcf = ParsedVcf()
                    oref = ref
                    parsed_vcf.ref = ref
                else:

                    #make some vars so easier to read
                    pos = split_line[1]
                    filter_flag = split_line[6]
                    var_call = split_line[4]
                    ref_call = split_line[3]

                    parsed_vcf.filter_flag[pos] = filter_flag

                    #split into ignore
                    if parsed_vcf.filter_flag[pos] != 'PASS':
                        parsed_vcf.bad_pos.append(pos)
                    else:
                        parsed_vcf.good_var.append(pos)

                    parsed_vcf.var[pos] = var_call
                    parsed_vcf.ref_base[pos] = ref_call

                    # get depth
                    matchObj = re.match(r'.*DP=(\d*)', line)
                    try:
                        parsed_vcf.depth[pos] = matchObj.group(1)
                    except AttributeError:
                        parsed_vcf.depth[pos] = 0
                    # get map quality
                    matchObj = re.match(r'.*MQ=(\d*\.\d*)', line)
                    try:
                        parsed_vcf.qual[pos] = matchObj.group(1)
                    except AttributeError:
                        parsed_vcf.qual[pos] = 0

                    #if a variant get some other things - including with it's a mix
                    if var_call != '.':
                        matchObvar = re.match(r'.*GT:AD:DP:GQ:PL\s+(.*)', line)
                        format_string = matchObvar.group(1).split(':')
                        parsed_vcf.hap_call[pos] = format_string[0]
                        parsed_vcf.hap_depth[pos] = format_string[2]
                        parsed_vcf.hap_qual[pos] = format_string[3]
                        ad_string = format_string[1].split(',')
                        parsed_vcf.hap_var_count[pos] = float(ad_string[1]) / float(parsed_vcf.depth[pos])
                        if parsed_vcf.hap_var_count[pos] < self.ad_cutoff:
                            parsed_vcf.mixed_positions.append(int(pos))


                 
            
        #add the last vcf
        self.parsed_vcf_container.append(parsed_vcf)
        #close file
        openfile.close()

        #calculate total number of mixed positions 
        self.number_mixed_positions = 0
        for p_vcf in self.parsed_vcf_container:
            self.number_mixed_positions += len(p_vcf.mixed_positions)

        #CHANGE we need Pheonix to calculate VCF depth
        self.depth_average = '31'        


    def define_class_variables_and_make_output_files(self, args):
        try:
            self.sample_name = os.path.basename(args.vcf[0]).split(os.extsep)[0]
        except AttributeError:
            self.sample_name = os.path.basename(args.fastqs[0]).split(os.extsep)[0]
        try:
            if os.path.exists(os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fa')):
                self.ref_genome_path = os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fa')
            elif os.path.exists(os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fasta')):
                self.ref_genome_path = (os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fasta'))
            elif os.path.exists(os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome)):
                self.ref_genome_path = (os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome))
        except IOError:
            sys.stderr.write('Cant find reference genome %s' % self.reference_genome)

        self.make_tmp_dir(args)
        #self.sorted_bamfile = os.path.join(self.tmp_dir, self.sample_name + '.sorted' + '.bam')
        self.vcf_filehandle = os.path.join(self.tmp_dir, os.path.pardir, '{0}.filtered.vcf'.format(self.sample_name))

    def check_reference_bwa_indexed(self):
        indices = ['amb', 'ann', 'bwt', 'pac', 'sa']
        for i in indices:
            if os.path.exists(self.ref_genome_path + '.' + i):
                pass
            else:
                os.system('bwa index %s' % self.ref_genome_path)
                # process = subprocess.Popen(['bwa', 'index', self.ref_genome_path],
                #                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for i in indices:
            if os.path.exists(self.ref_genome_path + '.' + i):
                pass
            else:
                sys.stderr.write('Reference not indexed, you need to index with `bwa index <ref_genome.fa>`')
                sys.exit()

    def run_phoenix(self,args):
        self.check_reference_bwa_indexed()
        self.check_reference_gatk_indexed()
        os.system('python run_snp_pipeline.py -r1 %s -r2 %s -r %s -o %s -m bwa -v gatk --sample-name %s --filters mq_score:%s,min_depth:%s,ad_ratio:%s' % args.fastqs[0], args.fastqs[1], self.reference_genome, self.tmp_dir, self.sample_name, self.mq_cutoff, self.depth_cutoff,self.ad_cutoff)



    def check_reference_gatk_indexed(self):
        indicies = ['dict', 'fa.fai']
        for i in indicies:
            if os.path.exists(os.path.splitext(self.ref_genome_path)[0] + '.' + i):
                pass
            else:
                picard_dict_path = os.path.splitext(self.ref_genome_path)[0]
                os.system('java -jar CreateSequenceDictionary.jar R= %s O= %s.dict'
                          % (self.ref_genome_path, picard_dict_path))
                os.system('samtools faidx %s' % self.ref_genome_path)
        for i in indicies:
            if os.path.exists(os.path.splitext(self.ref_genome_path)[0] + '.' + i):
                pass
            else:
                sys.stderr.write('Reference not indexed for GATK, you need to index with according to '
                                 'https://www.broadinstitute.org/gatk/guide/article?id=1601')
                sys.exit()


    def read_rec_file(self, rec_file):
        try:
            openfile = open(rec_file, 'r')
        except:
            print (rec_file + " not found ... ")
            sys.exit()
        rec_list = []
        for line in openfile:
            if line[0].isdigit():
                temp = (line.strip()).split('\t')
                rec_range = range((int(temp[0]) - 1), (int(temp[1]) - 1))
                rec_list = set(rec_list) | set(rec_range)
        self.rec_list = rec_list
