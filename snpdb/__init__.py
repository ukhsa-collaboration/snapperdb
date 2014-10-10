__author__ = 'flashton'

import logging
import psycopg2, psycopg2.extras
import sys
import os
from Bio import SeqIO
import datetime
import inspect
import SnapperDB
import errno
from gbru_vcf import Vcf
import pickle

class SNPdb:
    """
    Some really insteresting documnetation about the awesome class that this is!

    Examples
    --------
    To use this in your code do the following:

         # >>> from blah import SNPdb
         # >>> db = SNPdb(config="/path")
         # >>> db.get_matrix()
    """

    path_to_config = None
    """Path to the config"""
    snpdb_name = None
    """Name of the database"""

    reference_genome = None
    pg_uname = None
    pg_pword = None
    pg_host = None
    conn_string = None
    ## concious design decision to make cutoffs part of the SNPdb module rather than vcf, as should be consistent within a db
    depth_cutoff = None
    mq_cutoff = None
    ad_cutoff = None


    def __init__(self, config_dict):
        """
        Constructor for the SNPdb class. If config is specified, then it will always be used to initialise
            the connection.

        Parameters
        ----------
        config: str
            Path to the config file.
        uname: str
            Username to be used for the db.

        See also
        --------

        __N.B.__ Some interesting note goes here.

        Notes
        -----

        """

        self.parse_config_dict(config_dict)
        self.snpdb_conn = None

        self.ref_genome_dir = SnapperDB.__ref_genome_dir__

        if os.path.exists(self.ref_genome_dir):
            pass
        else:
            sys.stderr.write('Ref genome dir %s not found\n' % self.ref_genome_dir)
            sys.exit()



    def parse_config_dict(self, config_dict):
        ## we loop through thusly in case not all these things are in the config
        for attr in config_dict:
            if attr == 'reference_genome':
                self.reference_genome = config_dict[attr]
            if attr == 'depth_cutoff':
                self.depth_cutoff = config_dict[attr]
            if attr == 'mq_cutoff':
                self.mq_cutoff = config_dict[attr]
            if attr == 'ad_cutoff':
                self.ad_cutoff = config_dict[attr]
            if attr == 'snpdb_name':
                self.snpdb_name = config_dict[attr]
            if attr == 'pg_uname':
                self.pg_uname = config_dict[attr]
            if attr == 'pg_pword':
                self.pg_pword = config_dict[attr]
            if attr == 'pg_host':
                self.pg_host = config_dict[attr]

    def mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def make_tmp_dir(self, args):
        self.tmp_dir = os.path.join(os.path.dirname(args.vcf[0]), 'snpdb', 'tmp')
        self.mkdir_p(self.tmp_dir)

    def define_class_variables_and_make_output_files(self, args, vcf):
        ## need to handle either vcf or fastqs
        try:
            vcf.sample_name = os.path.basename(args.vcf[0]).split(os.extsep)[0]
        except AttributeError:
            vcf.sample_name = os.path.basename(args.fastqs[0]).split(os.extsep)[0]
        vcf.ref_genome_path = os.path.join(SnapperDB.__ref_genome_dir__, self.reference_genome + '.fa')
        vcf.make_tmp_dir(args)
        vcf.sorted_bamfile = os.path.join(vcf.tmp_dir, vcf.sample_name + '.sorted' + '.bam')
        vcf.vcf_filehandle = os.path.join(vcf.tmp_dir, os.path.pardir, '{0}.vcf'.format(vcf.sample_name))

    def _write_conn_string(self):
        self.conn_string = 'host=\'{0}\' dbname={1} user=\'{2}\' password=\'{3}\''.format(self.pg_host, self.snpdb_name,
                                                                                          self.pg_uname, self.pg_pword)

    def _check_if_snpdb_exists(self):
        try:
            psycopg2.connect(self.conn_string)
            return True
        except psycopg2.OperationalError:
            return False

    def make_snpdb(self):
        does_snpdb_exist = self._check_if_snpdb_exists()
        if does_snpdb_exist == True:
            sys.stderr.write('This snpdb already exists\n')
        else:
            sys.stdout.write('The SNPdb {0} does not exist - running sql to  make snpdb\n'.format(self.snpdb_name))
            make_db_conn_string = 'host=\'{0}\' dbname=postgres user=\'{1}\' password=\'{2}\''.format(self.pg_host,
                                                                                          self.pg_uname, self.pg_pword)

            conn = psycopg2.connect(make_db_conn_string)
            conn.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)
            cur = conn.cursor()
            cur.execute('CREATE DATABASE {0}'.format(self.snpdb_name))
            cur.close()

            conn_string = 'host=\'{0}\' dbname={1} user=\'{2}\' password=\'{3}\''.format(self.pg_host, self.snpdb_name,
                                                                                         self.pg_uname, self.pg_pword)
            conn = psycopg2.connect(conn_string)
            cur = conn.cursor()
            sql_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'clean_snpdb.backup')
            cur.execute(open(sql_script, 'r').read())
            conn.commit()
            conn.close()


    def check_len_vcf(self, vcf):
        vcf_len = len(vcf.depth)
        ref_genome_path = os.path.join(self.ref_genome_dir, self.reference_genome + '.fa')
        with open(ref_genome_path, 'r') as fi:
            ref_fasta = SeqIO.parse(fi, 'fasta')
            if len(ref_fasta) == vcf_len:
                pass
            else:
                sys.stderr.write('VCF length and reference fasta length are not the same\n')
                sys.exit()

    def check_duplicate(self, vcf):
        dup = False
        dict_cursor = self.snpdb_conn.cursor(cursor_factory = psycopg2.extras.DictCursor)
        dict_cursor.execute("select distinct(name) FROM strains_snps where name = \'%s\'" % vcf.sample_name)
        for row in dict_cursor:
            dup = True
        dict_cursor.close()
        return dup

    def add_depth_to_strain_stats(self, vcf):
        time_now = datetime.datetime.now()
        time_now = str(time_now)
        insert_statement = 'INSERT INTO strain_stats (name, av_cov, time_of_upload) VALUES (\'%s\', %s, \'%s\')' % (
            vcf.sample_name, vcf.depth_average, time_now)
        cur = self.snpdb_conn.cursor()
        cur.execute(insert_statement)
        self.snpdb_conn.commit()
        cur.close()
        pass

    def add_new_variants(self, pos, contig, good_var):
        var_base = good_var[pos]
        # fasta starts at 0, snpdb pos is worked out from a reference that starts at 1 (position in snpdb)
        ref_base = contig[int(pos) - 1]

        cursor = self.snpdb_conn.cursor()
        cursor.execute("insert into variants (pos, var_base, ref_base) VALUES (%s,\'%s\',\'%s\')" % (pos, var_base, ref_base))
        self.snpdb_conn.commit()

        cursor.execute('select currval(\'variants_id_seq\')')
        seq_id = cursor.fetchall()
        seq_id = seq_id[0][0]
        return seq_id

    def add_to_snpdb(self, vcf):
        ref_genome_path = os.path.join(self.ref_genome_dir, self.reference_genome + '.fa')
        if os.path.exists(ref_genome_path):
            pass
        else:
            sys.stderr.write('Could not find {0}, check your file extension (needs to be .fa)\n'.format(ref_genome_path))

        with open(ref_genome_path, 'r') as fi:
            ref_fasta = SeqIO.parse(fi, 'fasta')
            contig = ''
            i = 0
            for each in ref_fasta:
                i += 1
                if i != 1:
                    raise Exception
                elif i == 1:
                    contig = each.seq

            var_dic = {}
            dict_cursor = self.snpdb_conn.cursor(cursor_factory = psycopg2.extras.DictCursor)
            dict_cursor.execute('SELECT * FROM variants')
            for row in dict_cursor:
                #print row
                #print row['pos']
                if row['pos'] in var_dic:
                    var_dic[row['pos']][row['var_base']] = row['id']
                else:
                    var_dic[row['pos']] = {}
                    var_dic[row['pos']][row['var_base']] = row['id']

            var_db_list = []

            for pos in vcf.good_var:
                if int(pos) not in var_dic:
                    seq_id = self.add_new_variants(pos, contig, vcf.good_var)
                    var_db_list.append(seq_id)
                elif vcf.good_var[pos] not in var_dic[int(pos)]:
                    seq_id = self.add_new_variants(pos, contig, vcf.good_var)
                    var_db_list.append(seq_id)
                else:
                    var_db_list.append(var_dic[int(pos)][vcf.good_var[pos]])

            time_now = datetime.datetime.now()
            time_now = str(time_now)
            insert_statement = "insert into strains_snps (name, time_of_upload, variants_id, ignored_pos) VALUES (" \
                               "\'" + vcf.sample_name + "\', \'" + time_now + "\' ,\'{"

            for var_id in var_db_list:
                insert_statement += str(var_id) + ","
            insert_statement = insert_statement[:-1]
            insert_statement += "}\',\'{"

            for pos in vcf.bad_pos:
                insert_statement += str(pos) + ","
            insert_statement = insert_statement[:-1]

            insert_statement += "}\')"

            cursor = self.snpdb_conn.cursor()
            cursor.execute(insert_statement)
            self.snpdb_conn.commit()
            cursor.close()
            self.snpdb_conn.close()

    def snpdb_upload(self, vcf):
        if self.check_duplicate(vcf) == False:
            self.add_depth_to_strain_stats(vcf)
            self.add_to_snpdb(vcf)
        elif self.check_duplicate(vcf) == True:
            sys.stderr.write('Sample is already in SNPdb\n')


def vcf_to_db(args, config_dict, vcf):
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._write_conn_string()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)
    print snpdb.conn_string
    if inspect.stack()[0][3] == 'fastq_to_db':
        pass
    elif inspect.stack()[0][3] == 'vcf_to_db':
        ## there is no existing vcf class here, but there will definitely be a vcf, and there may be a pickle.
        vcf = Vcf()
        snpdb.define_class_variables_and_make_output_files(args, vcf)
        if os.path.exists(os.path.join(vcf.tmp_dir, vcf.sample_name  + '_bad_pos.pick')):
            bad_pos = pickle.load(open(os.path.join(vcf.tmp_dir, vcf.sample_name  + '_bad_pos.pick')))
            good_var = pickle.load(open(os.path.join(vcf.tmp_dir, vcf.sample_name  + '_good_var.pick')))
            vcf.good_var = good_var
            vcf.bad_pos = bad_pos
            snpdb.check_len_vcf(vcf)
            snpdb.snpdb_upload(vcf)
        else:
            vcf.parse_config_dict(config_dict)
            vcf.read_vcf()
            snpdb.snpdb_upload(vcf)

            '''
            Parse vcf and get good_var and ignored_pos
            '''

def make_snpdb(config_dict):
    snpdb = SNPdb(config_dict)
    snpdb._write_conn_string()
    snpdb.make_snpdb()

