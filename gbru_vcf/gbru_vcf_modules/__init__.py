__author__ = 'flashton'

import re
import os
import sys
import subprocess
try:
    import cPickle as pickle
except:
    import pickle


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

    def read_vcf(self, file_name, SNPdb):
        ## as the cutoffs are part of the SNPdb, need to pass the SNPdb instance
        try:
            openfile = open(file_name, 'r')
        except:
            print file_name + " not found ... "
            sys.exit()
        self.sample_name = os.path.splitext(os.path.basename(file_name))[0]

        #get query from file name
        temp = file_name.split('/')
        temp2 = temp[len(temp)-1].split('.')
        name = temp2[0]
        self.query = temp2[0]

        for line in openfile:
            if not re.match(r'^#',line):
                temp = line.split()
                #get reference
                self.ref = temp[0]
                #get pos
                pos = temp[1]
                #get depth
                matchObj  = re.match(r'.*DP=(\d*)',line)
                try:
                    self.depth[pos] = matchObj.group(1)
                except:
                    self.depth[pos] = 0
                #get map quality
                matchObj  = re.match(r'.*MQ=(\d*\.\d*)',line)
                try:
                    self.qual[pos] = matchObj.group(1)
                except:
                    self.qual[pos]  = 0
                #get var call
                var_call = temp[4]
                ref_call = temp[3]
                #if not wild type
                if var_call != '.':
                    matchObvar  = re.match(r'.*GT:AD:DP:GQ:PL\s+(.*)',line)
                    temp = matchObvar.group(1).split(':')
                    self.hap_call[pos] = temp[0]
                    self.hap_depth [pos] = temp[2]
                    self.hap_qual[pos] = temp[3]
                    temp2 = temp[1].split(',')
                    self.hap_var_count[pos] = float(temp2[1]) / float(self.depth[pos])
                    self.var[pos] = var_call
                    self.ref_base[pos] = ref_call

        self.bad_depth = self.return_positions_with_depth_less_than_cutoff(SNPdb.depth_cutoff)
        self.bad_qual = self.return_positions_with_mq_less_than_cutoff(SNPdb.mq_cutoff)
        self.bad_var, self.good_var = self.return_var_positions_with_depth_qual_less_and_above_than_cutoff(SNPdb.depth_cutoff,
                                                                                                           SNPdb.mq_cutoff,
                                                                                                           SNPdb.ad_cutoff)
        self.bad_pos = set(self.bad_depth) | set(self.bad_qual) | set(self.bad_var)

        self.get_average_and_sd_depth()

    def get_length_of_ref(self):
        ref_len = len(self.depth)
        return ref_len

    def get_average_and_sd_depth(self):
        ref_len = len(self.depth)
        total = 0
        for pos in self.depth:
            total = total + float(self.depth[pos])
        self.depth_average = float(total) / float(ref_len)

        sd_tot=0
        for pos in sorted(self.depth, key=int):
            sd_tot = sd_tot + (float(self.depth[pos]) - float(self.depth_average))**2

        self.depth_sd = (float(sd_tot)/float(ref_len))**0.5

    def return_positions_with_depth_less_than_cutoff(self,cutoff):
        bad_list = []
        for pos in self.depth:
            if float(self.depth[pos]) < float(cutoff):
                bad_list.append(pos)
        return bad_list

    def return_positions_with_mq_less_than_cutoff(self,cutoff):
        bad_list = []
        for pos in self.qual:
            if float(self.qual[pos]) < float(cutoff):
                bad_list.append(pos)
        return bad_list

    def return_var_positions_with_depth_qual_less_and_above_than_cutoff(self, depth_co, qual_co,var_count_co):
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

    def pickle_variants_and_ignored_pos(self, opts):
        opts.bad_pos_pick = os.path.join(os.path.split(opts.vcf)[0], 'tmp', '{0}_bad_pos.pick'.format(os.path.split(opts.vcf)[
            1]))
        opts.good_var_pick = os.path.join(os.path.split(opts.vcf)[0], 'tmp', '{0}_good_var.pick'.format(os.path.split(opts.vcf)[
            1]))
        with open(opts.bad_pos_pick, 'wb') as fo:
            pickle.dump(self.bad_pos, fo, -1)
        with open(opts.good_var_pick, 'wb') as fo:
            pickle.dump(self.good_var, fo, -1)

    def parse_config(self):
        try:
            with open(self.path_to_config, 'r') as fi:
                for line in fi.readlines():
                    if line.startswith('reference_genome'):
                        self.reference_genome = line.strip().split()[-1]

        except IOError:
            print 'Cannot find {0}'.format(self.path_to_config)
            sys.exit()

        self.ref_genome_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname( __file__ ))),
                                                                     'reference_genomes'))
        if os.path.exists(self.ref_genome_dir):
            pass
        else:
            sys.stderr.write('Ref genome dir not found\n')
            sys.exit()

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
        process = subprocess.Popen(['/Users/flashton/Programs/bwa-0.7.5a/bwa', 'mem', '-R', header, self.ref_genome_path,
                                    fastq_1, fastq_2],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write('Problem with bwa mapping\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()
        samfile = os.path.join(out_dir, self.sample_name + '.sam')
        with open(samfile, 'w') as fo:
            fo.write(stdout)

    def convert_sort_index(self, opts):
        samfile = os.path.join(opts.out_dir, self.sample_name + '.sam')
        bamfile = os.path.join(opts.out_dir, self.sample_name + '.bam')
        sorted_bamfile = os.path.join(opts.out_dir, self.sample_name + '.sorted')
        process = subprocess.Popen(['samtools', 'view', '-bS', '-o', bamfile, samfile], stderr = subprocess.PIPE)
        stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write('Problem with <convert>_sort_index\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()

        process = subprocess.Popen(['samtools', 'sort', bamfile, sorted_bamfile], stderr = subprocess.PIPE)
        stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write('Problem with convert_<sort>_index\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()

        sorted_bamfile = os.path.join(opts.out_dir, self.sample_name + '.sorted' + '.bam')
        process = subprocess.Popen(['samtools', 'index', sorted_bamfile], stderr = subprocess.PIPE)
        stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write('Problem with convert_sort_<index>\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()

    def make_sorted_bam(self, opts):
        self.parse_config()
        self.ref_genome_path = os.path.join(self.ref_genome_dir, self.reference_genome + '.fa')
        self.check_reference_bwa_indexed()
        self.run_bwa(opts.fastq_1, opts.fastq_2, opts.out_dir)
        self.convert_sort_index(opts)

    def check_reference_gatk_indexed(self):
        indicies = ['dict', 'fa.fai']
        for i in indicies:
            if os.path.exists(os.path.splitext(self.ref_genome_path)[0] + '.' + i):
                pass
            else:
                sys.stderr.write('Reference not indexed for GATK, you need to index with according to '
                                 'https://www.broadinstitute.org/gatk/guide/article?id=1601')
                sys.exit()

    def run_gatk(self, opts):
        sorted_bamfile = os.path.join(opts.out_dir, self.sample_name + '.sorted' + '.bam')
        vcf_file = os.path.join(opts.out_dir, os.path.pardir, '{0}.vcf'.format(self.sample_name))
        process = subprocess.Popen(['java', '-Xmx30g', '-jar',
                                    '/Users/flashton/Programs/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar', '-T',
                                    'UnifiedGenotyper', '-nt', '1', '-out_mode', 'EMIT_ALL_SITES', '-R', self.ref_genome_path,
                                    '-I', sorted_bamfile, '-o', vcf_file], stderr = subprocess.PIPE )
        stderr = process.communicate()

        if process.returncode != 0:
            sys.stderr.write('Problem with run_gatk\n')
            sys.stderr.write('{0}\n'.format(stderr))
            sys.exit()

    def make_vcf(self, opts):
        self.check_reference_gatk_indexed()
        self.run_gatk(opts)