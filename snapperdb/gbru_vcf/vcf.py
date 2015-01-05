__author__ = 'flashton'

import errno
import os
import pickle
import re
import subprocess
import sys

import snapperdb


class Vcf:
    def __init__(self):
        self.depth = {}
        self.qual = {}
        self.hap_qual = {}
        self.hap_depth = {}
        self.hap_call = {}
        self.hap_var_count = {}
        self.var = {}
        self.ref_base = {}
        self._query = str
        self._ref = str
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
        self.sorted_bamfile = None
        self.vcf_filehandle = None
        self.depth_cutoff = None
        self.mq_cutoff = None
        self.ad_cutoff = None
        self.number_mixed_positions = 0
        self.mixed_positions = []
        self.rec_list = []

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

    def read_vcf(self):

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
        self.sample_name = os.path.splitext(os.path.basename(self.vcf_filehandle))[0]

        # # bit pointless but maintains 'query' syntax of the original
        self.query = self.sample_name
        print self.query

        for line in openfile:
            if not re.match(r'^#', line):
                split_line = line.split()
                # get reference
                self.ref = split_line[0]
                # get pos
                pos = split_line[1]
                # get depth
                matchObj = re.match(r'.*DP=(\d*)', line)
                try:
                    self.depth[pos] = matchObj.group(1)
                except:
                    self.depth[pos] = 0
                # get map quality
                matchObj = re.match(r'.*MQ=(\d*\.\d*)', line)
                try:
                    self.qual[pos] = matchObj.group(1)
                except:
                    self.qual[pos] = 0
                # get var call
                var_call = split_line[4]
                ref_call = split_line[3]
                # if not wild type
                if var_call != '.':
                    matchObvar = re.match(r'.*GT:AD:DP:GQ:PL\s+(.*)', line)
                    format_string = matchObvar.group(1).split(':')
                    self.hap_call[pos] = format_string[0]
                    self.hap_depth [pos] = format_string[2]
                    self.hap_qual[pos] = format_string[3]
                    ad_string = format_string[1].split(',')
                    self.hap_var_count[pos] = float(ad_string[1]) / float(self.depth[pos])
                    if self.hap_var_count[pos] < self.ad_cutoff:
                        self.mixed_positions.append(pos)
                        # if pos not in self.rec_list:
                        #    self.number_mixed_positions += 1
                    self.var[pos] = var_call
                    self.ref_base[pos] = ref_call


        self.bad_depth = self.return_positions_with_low_depth(self.depth_cutoff)
        self.bad_qual = self.return_positions_with_low_mq(self.mq_cutoff)
        self.bad_var, self.good_var = self.return_bad_pos_good_vars(self.depth_cutoff,
                                                                                                           self.mq_cutoff,
                                                                                                           self.ad_cutoff)
        self.bad_pos = set(self.bad_depth) | set(self.bad_qual) | set(self.bad_var)

        self.get_average_and_sd_depth()
        openfile.close()

    def get_length_of_ref(self):
        ref_len = len(self.depth)
        return ref_len

    def get_average_and_sd_depth(self):
        ref_len = len(self.depth)
        total = 0
        for pos in self.depth:
            total += float(self.depth[pos])
        self.depth_average = float(total) / float(ref_len)

        sd_tot = 0
        for pos in sorted(self.depth, key=int):
            sd_tot = sd_tot + (float(self.depth[pos]) - float(self.depth_average)) ** 2

        self.depth_sd = (float(sd_tot) / float(ref_len)) ** 0.5

    def return_positions_with_low_depth(self, cutoff):
        bad_list = []
        for pos in self.depth:
            if float(self.depth[pos]) < float(cutoff):
                bad_list.append(pos)
        return bad_list

    def return_positions_with_low_mq(self, cutoff):
        bad_list = []
        for pos in self.qual:
            if float(self.qual[pos]) < float(cutoff):
                bad_list.append(pos)
        return bad_list

    def return_bad_pos_good_vars(self, depth_co, qual_co, var_count_co):
        bad_list = []
        good_dict = {}
        for pos in self.hap_qual:
            if float(self.hap_qual[pos]) < float(qual_co):
                bad_list.append(pos)
            elif float(self.hap_depth[pos]) < float(depth_co):
                bad_list.append(pos)
            elif self.hap_call[pos] != '1/1':
                bad_list.append(pos)
            elif float(self.hap_var_count[pos]) < float(var_count_co):
                bad_list.append(pos)
            else:
                good_dict[pos] = self.var[pos]
        return bad_list, good_dict

    def pickle_variants_and_ignored_pos(self, args):
        args.bad_pos_pick = os.path.join(self.tmp_dir, '{0}_bad_pos.pick'.format(os.path.split(self.sample_name)[
            1]))
        args.good_var_pick = os.path.join(self.tmp_dir, '{0}_good_var.pick'.format(os.path.split(self.sample_name)[
            1]))
        with open(args.bad_pos_pick, 'wb') as fo:
            pickle.dump(self.bad_pos, fo, -1)
        with open(args.good_var_pick, 'wb') as fo:
            pickle.dump(self.good_var, fo, -1)

    def define_class_variables_and_make_output_files(self, args):
        try:
            self.sample_name = os.path.basename(args.vcf[0]).split(os.extsep)[0]
        except AttributeError:
            self.sample_name = os.path.basename(args.fastqs[0]).split(os.extsep)[0]

        self.ref_genome_path = os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fa')
        self.make_tmp_dir(args)
        self.sorted_bamfile = os.path.join(self.tmp_dir, self.sample_name + '.sorted' + '.bam')
        self.vcf_filehandle = os.path.join(self.tmp_dir, os.path.pardir, '{0}.vcf'.format(self.sample_name))

    def check_reference_bwa_indexed(self):
        indices = ['amb', 'ann', 'bwt', 'pac', 'sa']
        for i in indices:
            if os.path.exists(self.ref_genome_path + '.' + i):
                pass
            else:
                process = subprocess.Popen(['bwa', 'index', self.ref_genome_path],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for i in indices:
            if os.path.exists(self.ref_genome_path + '.' + i):
                pass
            else:
                sys.stderr.write('Reference not indexed, you need to index with `bwa index <ref_genome.fa>`')
                sys.exit()

    def run_bwa(self, fastq_1, fastq_2, out_dir):
        header = '@RG\tID:1\tSM:%s' % self.sample_name
        process = subprocess.Popen(['bwa', 'mem', '-R', header, self.ref_genome_path,
                                    fastq_1, fastq_2],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write('Problem with bwa mapping\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()
        samfile = os.path.join(out_dir, '{0}.sam'.format(self.sample_name))
        with open(samfile, 'w') as fo:
            fo.write(stdout)

    def convert_sort_index(self, args):
        samfile = os.path.join(self.tmp_dir, self.sample_name + '.sam')
        bamfile = os.path.join(self.tmp_dir, self.sample_name + '.bam')
        # # have to make tmp_sorted_bamfile because samtools sort takes filename without .bam as output handle
        tmp_sorted_bamfile = os.path.join(self.tmp_dir, self.sample_name + '.sorted')
        process = subprocess.Popen(['samtools', 'view', '-bS', '-o', bamfile, samfile], stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write('Problem with <convert>_sort_index\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()

        process = subprocess.Popen(['samtools', 'sort', bamfile, tmp_sorted_bamfile], stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write('Problem with convert_<sort>_index\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()

        process = subprocess.Popen(['samtools', 'index', self.sorted_bamfile], stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write('Problem with convert_sort_<index>\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()

    def make_sorted_bam(self, args):
        self.check_reference_bwa_indexed()
        self.run_bwa(args.fastqs[0], args.fastqs[1], self.tmp_dir)
        self.convert_sort_index(args)

    def check_reference_gatk_indexed(self):
        indicies = ['dict', 'fa.fai']
        for i in indicies:
            if os.path.exists(os.path.splitext(self.ref_genome_path)[0] + '.' + i):
                pass
            else:
                sys.stderr.write('Reference not indexed for GATK, you need to index with according to '
                                 'https://www.broadinstitute.org/gatk/guide/article?id=1601')
                sys.exit()

    def run_gatk(self, args):
        process = subprocess.Popen(['java', '-Xmx30g', '-jar',
                                    '/phengs/hpc_software/gatk/2.6.5/GenomeAnalysisTK.jar', '-T',
                                    'UnifiedGenotyper', '-nt', '1', '-out_mode', 'EMIT_ALL_SITES', '-R', self.ref_genome_path,
                                    '-I', self.sorted_bamfile, '-o', self.vcf_filehandle], stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write('Problem with run_gatk\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()

    def make_vcf(self, args):
        self.check_reference_gatk_indexed()
        self.run_gatk(args)

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