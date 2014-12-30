__author__ = 'flashton'

from datetime import datetime
import errno
import glob
import inspect
import os
import pickle
import re
import sys
import logging
from Bio import SeqIO
import psycopg2, psycopg2.extras

import snapperdb
from snapperdb.gbru_vcf import Vcf


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
    # # concious design decision to make cutoffs part of the SNPdb module rather than vcf, as should be consistent within a db
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
        self.strains_snps = {}
        self.parse_config_dict(config_dict)
        self.snpdb_conn = None
        self.ref_genome_dir = snapperdb.__ref_genome_dir__
        if os.path.exists(self.ref_genome_dir):
            pass
        else:
            sys.stderr.write('Ref genome dir %s not found\n' % self.ref_genome_dir)
            sys.exit()

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
        # # need to handle either vcf or fastqs
        try:
            vcf.sample_name = os.path.basename(args.vcf[0]).split(os.extsep)[0]
        except AttributeError:
            vcf.sample_name = os.path.basename(args.fastqs[0]).split(os.extsep)[0]
        vcf.ref_genome_path = os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fa')
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
        print ref_genome_path
        fi = open(ref_genome_path)
        ref_fasta = SeqIO.read(fi, 'fasta')
        print len(ref_fasta.seq), len(vcf.depth)
        if len(ref_fasta.seq) == vcf_len:
            pass
        else:
            sys.stderr.write('VCF length and reference fasta length are not the same\n')
            sys.exit()
        fi.close()

    def check_duplicate(self, vcf):
        dup = False
        dict_cursor = self.snpdb_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        dict_cursor.execute("select distinct(name) FROM strains_snps where name = \'%s\'" % vcf.sample_name)
        for row in dict_cursor:
            dup = True
        dict_cursor.close()
        return dup

    def add_depth_to_strain_stats(self, vcf):
        time_now = datetime.now()
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
            dict_cursor = self.snpdb_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            dict_cursor.execute('SELECT * FROM variants')
            for row in dict_cursor:
                # print row
                # print row['pos']
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

            time_now = datetime.now()
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

    # # functions below here are for querying the snpdb

    def get_background(self, cur, strain_list, args):
        sql = "select distinct(" + args.back_flag + ") from strain_clusters"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            sql2 = "select sc.name from strain_clusters sc , strain_meta sm where sm.name = sc.name and " + args.back_flag + "=" + str(
                row[0]) + " "
            if args.meta_flag != 'N':
                temp = args.meta_flag.split(',')
                for meta in temp:
                    temp2 = meta.split(':')
                    sql2 = sql2 + "and " + temp2[0] + "=\'" + str(temp2[1]) + "\' "
                sql2 = sql2 + " limit 1"
            # print sql2
            cur.execute(sql2)
            rows2 = cur.fetchone()
            if rows2:
                strain_list.append(rows2[0])
        return strain_list

    def add_strains_to_sql_co(self, sql, strain_list, co, name):
        for strain in strain_list:
            sql += "name =\'" + strain + "\' or "
        sql = sql[:-4]
        sql += ") and icount(" + name + ") < " + co
        return sql

    def get_all_good_ids(self, cur, strain_list, snp_co):
        strain_snps = {}
        totlist = []
        # # this could be replaced by where like any (array[strain_list]) - I think
        sql = " select variants_id, name, icount(variants_id) as count from strains_snps where ("
        sql = self.add_strains_to_sql_co(sql, strain_list, snp_co, "variants_id")
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            totlist = set(totlist) | set(row[0])
            strain_snps[row[1]] = row[0]
        return totlist, strain_snps

    def get_variants(self, cur):
        variant_container = {}
        pos_2_id_list = {}
        sql = "select pos , id, ref_base, var_base, gene,product, amino_acid from variants"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            variant = Variant()
            variant.pos = row[0]
            variant.id = row[1]
            variant.ref_base = row[2]
            variant.var_base = row[3]
            variant.gene = row[4]
            variant.product = row[5]
            variant.amino_acid = row[6]
            variant_container[row[1]] = variant
            if row[1] in self.goodids:
                if row[0] in pos_2_id_list:
                    pos_2_id_list[row[0]].append(row[1])
                else:
                    pos_2_id_list[row[0]] = []
                    pos_2_id_list[row[0]].append(row[1])
        return variant_container, pos_2_id_list

    def get_bad_pos_for_strain_get_the_vars(self, cur, strain, totlist):
        strain_ig = {}
        sql = "select ignored_pos, icount(ignored_pos), name as count from strains_snps where name =\'" + strain + "\'"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            totlist = set(totlist) | set(row[0])
            strain_ig = row[0]
        return totlist, strain_ig

    def make_consensus(self, ref_seq, ref_flag, cur):
        fasta = {}
        var_id_list = []
        var_look = {}
        n_look = {}
        badlist = []
        for strain in self.strains_snps:
            badlist, strains_ig = self.get_bad_pos_for_strain_get_the_vars(cur, strain, badlist)
            if strain not in fasta:
                fasta[strain] = ref_seq[:]
            for ids in self.strains_snps[strain]:
                if self.variants[ids].var_base != fasta[strain][self.variants[ids].pos - 1]:
                    fasta[strain][self.variants[ids].pos - 1] = self.variants[ids].var_base
                    if self.variants[ids].pos:
                        if self.variants[ids].pos in var_look:
                            var_look[self.variants[ids].pos] = var_look[self.variants[ids].pos] + 1
                            if ids not in var_id_list:
                                var_id_list.append(ids)
                        else:
                            var_look[self.variants[ids].pos] = 1
                            var_id_list.append(ids)
            for bad_ids in strains_ig:
                fasta[strain][bad_ids - 1] = 'N'
                n_look[bad_ids] = 'N'
        if ref_flag == 'Y':
            fasta[self.reference_genome] = ref_seq
        return fasta, var_look, n_look, badlist, var_id_list

    def calc_matrix(self):
        diff_matrix = {}
        for strain1 in self.fasta:
            diff_matrix[strain1] = {}
            for strain2 in self.fasta:
                if strain2 in diff_matrix and strain1 not in diff_matrix[strain2]:
                    diff_matrix[strain1][strain2] = 0
                    for var in self.var_look:
                        if self.fasta[strain1][var - 1] != self.fasta[strain2][var - 1] and self.fasta[strain1][var - 1] != 'N' and self.fasta[strain2][var - 1] != 'N':                                 diff_matrix[strain1][strain2] += 1

        return diff_matrix

    def parse_args_for_get_the_snps(self, args, cur, strain_list, ref_seq):
        # # this populates class variables specific to querying the SNPdb
        if args.back_flag != 'N':
            print "###  Getting background strains:" + str(datetime.now())
            strain_list = self.get_background(cur, strain_list, args)

        print "###  Getting good positions:" + str(datetime.now())
        self.goodids, self.strains_snps = self.get_all_good_ids(cur, strain_list, args.snp_co)

        print "Variable positions: " + str(len(self.goodids))

        if len(self.goodids) == 0:
            print "###  No variable positions found: EXITING " + str(datetime.now())
            sys.exit()
        print str(len(self.strains_snps)) + " strains used out of " + str(len(strain_list))
        print "###  Getting variants:" + str(datetime.now())
        self.variants, self.posIDMap = self.get_variants(cur)
        print "###  Making consensus:" + str(datetime.now())
        self.fasta, self.var_look, self.n_look, self.badlist, self.var_id_list = self.make_consensus(ref_seq, args.ref_flag, cur)
        print "Ignored positions: " + str(len(self.badlist))
        if args.mat_flag == 'Y':
            print "###  Making Distance Matrix: " + str(datetime.now())
            self.matrix = self.calc_matrix()

    def print_fasta(self, out, flag, ref_flag, rec_list):
        f = open(out + '.fa', 'w')
        for strain in self.fasta:
            f.write(">" + strain + "\n")
            if flag == 'W':
                for i, seq in enumerate(self.fasta[strain]):
                    f.write(seq)
            elif flag == 'A':
                for i in sorted(self.var_look):
                    if i not in rec_list:
                        if ref_flag == 'Y':
                            f.write(self.fasta[strain][i - 1])
                        elif self.var_look[i] != len(self.strains_snps):
                            f.write(self.fasta[strain][i - 1])
            elif flag == 'C':
                for i in sorted(self.var_look):
                    if i not in self.n_look:
                        if i not in rec_list:
                            if ref_flag == 'Y':
                                f.write(self.fasta[strain][i - 1])
                            elif self.var_look[i] != len(self.strains_snps):
                                f.write(self.fasta[strain][i - 1])
            f.write("\n")

    def print_matrix(self, out):
        f = open(out + '.matrix', 'w')
        for strain1 in self.matrix:
            for strain2 in self.matrix[strain1]:
                if strain1 != strain2:
                    f.write(strain1 + "\t" + strain2 + "\t" + str(self.matrix[strain1][strain2]) + "\n")

    def print_vars(self, out, flag, rec_list, ref_flag):
        f = open(out + '.variants', 'w')
        for pos in sorted(self.var_look):
            var_ids = self.posIDMap[pos]
            for var_id in var_ids:
                if var_id in self.var_id_list:
                    if flag == 'W':
                        f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].gene) + "\n")
                    elif flag == 'A':
                        if pos not in rec_list:
                            if ref_flag == 'Y':
                                f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                            elif self.var_look[pos] != len(self.strains_snps):
                                f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")

                    elif flag == 'C':
                        if pos not in self.n_look:
                            if pos not in rec_list:
                                if ref_flag == 'Y':
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                                elif self.var_look[pos] != len(self.strains_snps):
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")

    # # functions below here are from update_distance_matrix

    def get_strains(self, cur):
        strain_list = []
        sql = " select name from strain_stats where ignore is NULL"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            strain_list.append(row[0])

        update_strain = []
        for strain1 in strain_list:
            sql = "select * from dist_matrix where strain1 = '%s' or strain2 = '%s' limit 1" % (strain1, strain1)
            cur.execute(sql)
            row = cur.fetchall()
            if not row:
                update_strain.append(strain1)

        return strain_list, update_strain

    def parse_args_for_update_matrix(self, cur, snp_co, strain_list):
        # # this populates class variables specific to querying the SNPdb

        print "###  Getting good positions:" + str(datetime.now())
        self.goodids, self.strains_snps = self.get_all_good_ids(cur, strain_list, snp_co)
        print "###  Getting variants:" + str(datetime.now())
        self.variants, self.posIDMap = self.get_variants(cur)

    def get_bad_pos_for_strain_update_matrix(self, cur, strain):
        strain_ig = []
        sql = "select ignored_pos, icount(ignored_pos), name from strains_snps where name =\'" + strain + "\'"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            strain_ig = row[0]
        return strain_ig

    def check_matrix(self, cur, data_list, update_strain):
        seen_strain = []
        for strain1 in update_strain:
            seen_strain.append(strain1)
            print "Populating matrix for: " + strain1
            strain1_good_var = self.strains_snps[strain1]
            strain1_ig_pos = self.get_bad_pos_for_strain_update_matrix(cur, strain1)
            for strain2 in data_list:
                if strain1 != strain2 and strain2 not in seen_strain:
                    strain2_good_var = self.strains_snps[strain2]
                    strain2_ig_pos = self.get_bad_pos_for_strain_update_matrix(cur, strain2)
                    # getunion of bad_pos
                    all_bad_pos = set(strain1_ig_pos) | set(strain2_ig_pos)

                    # getsymmetric difference of variants
                    all_var = set(strain1_good_var) ^ set(strain2_good_var)
                    diff = 0
                    for var_id in all_var:
                        if self.variants[var_id].pos not in all_bad_pos:
                            diff = diff + 1
                    # add to db
                    sql2 = "insert into dist_matrix (strain1, strain2, snp_dist) VALUES (\'%s\',\'%s\',%s)" % (strain1, strain2, diff)
                    cur.execute(sql2)
                    self.snpdb_conn.commit()

    def chunks(self, l, n):
        for i in xrange(0, len(l), n):
            yield l[i:i + n]

    def write_qsubs_to_check_matrix(self, args, strain_list, short_strain_list, update_strain, snpdb):
        '''
        how to call this particular function as a qsub - need access to check matrix on command line.
        1. write lists
        2. for each in lists, qsub
        '''
        home_dir = os.path.expanduser('~')
        logs_dir = os.path.join(home_dir, 'logs')
        this_dir = os.path.dirname(os.path.realpath(__file__))
        snapperdb_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

        for i, each in enumerate(self.chunks(update_strain, args.hpc)):
            with open('{0}/update_list_{1}'.format(this_dir, i), 'w') as fo:
                for x in each:
                    fo.write(x + '\n')
        with open('{0}/short_strain_list'.format(this_dir), 'w') as fo:
            for x in short_strain_list:
                fo.write(x + '\n')

        with open('{0}/strain_list'.format(this_dir), 'w') as fo:
            for x in strain_list:
                fo.write(x + '\n')

        res = sorted(glob.glob('{0}/update_list*'.format(this_dir)))

        for i, update_list in enumerate(res):
            command = ('#! /bin/bash\n'
                       '#$ -o {0}/check_matrix.stdout\n'
                       '#$ -e {0}/check_matrix.stderr\n'
                       '#$ -m e\n'
                       '#$ -wd {1}\n'
                       '#$ -N up_mat_{2}_{3}\n\n'
                       '. /etc/profile.d/modules.sh\n'
                       'module load {7}/.module_files/snapperdb/1-0\n'
                       'python SnapperDB_main.py'
                       ' qsub_to_check_matrix -c {4}'
                       ' -l {5}/strain_list'
                       ' -s {5}/short_strain_list'
                       ' -u {6}\n'.format(logs_dir, snapperdb_dir, snpdb, i, args.config_file, this_dir, update_list, home_dir))

            with open('{0}/update_matrix_{1}.sh'.format(this_dir, i), 'w') as fo:
                fo.write(command)
            os.system('qsub {0}/update_matrix_{1}.sh'.format(this_dir, i))


            # os.system('chmod u+x {0}/update_matrix_{1}.sh'.format(this_dir, i))
            # os.system('{0}/update_matrix_{1}.sh'.format(this_dir, i))

    def sweep_matrix(self):
        '''
        if the total number of strains in snpdb is n, check that for each strain there is n(n+1)/2 (or is ithis n-1?)
        entries in matrix. if there isn't, run some version of check matrix.

        or

        for each strain in the original update_strain list, fill in the matrix, where seen strain includes everything in the
        database that isn't in update strains.

        1. wait until all qsub_check_matrix jobs finished
        '''

        pass

    # ## all functions below here are to do with the clustering

    def get_input(self, cur):
        profile_dict = {}
        sql = "select strain1, strain2, snp_dist from dist_matrix"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            if row[0] not in profile_dict:
                profile_dict[row[0]] = {}
                profile_dict[row[0]][row[1]] = row[2]
            else:
                profile_dict[row[0]][row[1]] = row[2]
            if row[1] not in profile_dict:
                profile_dict[row[1]] = {}
                profile_dict[row[1]][row[0]] = row[2]
            else:
                profile_dict[row[1]][row[0]] = row[2]
        return profile_dict

    def get_cutoffs(self, cur):
        co = []
        cluster_dict = {}
        cluster_strain_list = {}
        # get columns names
        sql = "SELECT * FROM information_schema.columns WHERE table_schema = 'public' AND table_name = 'strain_clusters'"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            if row[3][0] == "t":
                co.append(row[3][1:])
        co = sorted(co, key=int)
        return co

    def make_links(self, profile_dict, co):
        clusters = {}
        for i, strain1 in enumerate(profile_dict):
            clusters[i] = []
            clusters[i].append(strain1)
            for strain2 in profile_dict[strain1]:
                if int(profile_dict[strain1][strain2]) <= int(co):
                    clusters[i].append(strain2)
        return clusters

    def make_links_update(self, profile_dict, co, cluster_dict):
        clusters = {}
        seen_strain = []
        for i, cluster in enumerate(sorted(cluster_dict, key=int)):
            clusters[cluster] = cluster_dict[cluster]
            for strain1 in cluster_dict[cluster]:
                seen_strain.append(strain1)
                for strain2 in profile_dict[strain1]:
                    if int(profile_dict[strain1][strain2]) <= int(co) and strain2 not in cluster_dict[cluster]:
                        clusters[cluster].append(strain2)
                        seen_strain.append(strain2)

        cluster_count = int(cluster)
        for i, strain1 in enumerate(profile_dict):
                if strain1 not in seen_strain:
                    cluster_count = cluster_count + 1
                    clusters[cluster_count] = []
                    clusters[cluster_count].append(strain1)
                    for strain2 in profile_dict[strain1]:
                        if int(profile_dict[strain1][strain2]) <= int(co) and strain2 not in seen_strain:
                            clusters[cluster_count].append(strain2)

        return clusters

    def define_clusters(self, slvs):
        clusters = {}
        for each in slvs:
            clusters[each] = set(slvs[each])
            for each2 in slvs:
                    int_set = set(slvs[each]) & set(slvs[each2])
                    if len(int_set) >= 1:
                        clusters[each] = set(clusters[each]) | set(slvs[each2])
                        # can we merge with any other clusters
                        for made_clusters in clusters:
                            int_set = set(clusters[each]) & set(clusters[made_clusters])
                            if len(int_set) >= 1:
                                clusters[each] = set(clusters[made_clusters]) | set(clusters[each])
                                clusters[made_clusters] = clusters[each]
        return clusters

    def remove_duplicate_clusters(self, clusters):
        clean_clusters = {}
        for cluster1 in clusters:
            clean_clusters[tuple(sorted(clusters[cluster1]))] = 1
        return clean_clusters

    def remove_duplicate_clusters_update(self, clusters):
        clean_clusters = {}
        for cluster1 in clusters:
            if tuple(sorted(clusters[cluster1])) not in clean_clusters:
                clean_clusters[tuple(sorted(clusters[cluster1]))] = cluster1
            elif int(cluster1) < int(clean_clusters[tuple(sorted(clusters[cluster1]))]):
                clean_clusters[tuple(sorted(clusters[cluster1]))] = cluster1
        return clean_clusters

    def print_slv_clusters(self, clusters, levels):
        strain_list = {}
        print "#\t" + str(levels)
        for co in (sorted(clusters, key=clusters.get)):
            # print co
            for i, cluster in enumerate(clusters[co]):
                for strain in list(cluster):
                    if strain not in strain_list:
                        strain_list[strain] = []
                    strain_list[strain].append(i + 1)
        # print strain_list
        for strain in (sorted(tuple(strain_list), key=strain_list.get)):
            hier = ""
            for clust in strain_list[strain]:
                hier = hier + str(clust) + "."
            print strain + "\t",
            print hier[:-1]

    def add_clusters_to_table(self, cur, clusters, levels):
        strain_list = {}
        for i, l in enumerate(levels):
            levels[i] = 't%s' % l
        print "#\t" + str(levels)
        for co in (sorted(clusters, key=clusters.get)):
            # print co
            for i, cluster in enumerate(clusters[co]):
                for strain in list(cluster):
                    if strain not in strain_list:
                        strain_list[strain] = []
                    strain_list[strain].append(i + 1)


        for strain in (sorted(tuple(strain_list), key=strain_list.get)):
            sql = 'INSERT INTO strain_clusters (name'
            for l in levels:
                sql = sql + ', %s' % l

            sql = sql + ') VALUES (\'%s\'' % strain

            for i in strain_list[strain]:
                sql = sql + ', %s' % i
            # sql = sql[:-2]
            sql = sql + ')'
            cur.execute(sql)
            self.snpdb_conn.commit()

    def add_clusters_to_existing_table(self, clusters, profile_dict, levels, cur, cluster_strain_list):
        strain_list = {}
        for co in (sorted(clusters, key=clusters.get)):
            for i, cluster in enumerate(clusters[co]):
                for strain in list(cluster):
                    if strain not in strain_list:
                        strain_list[strain] = []
                    strain_list[strain].append(int(clusters[co][cluster]))
        for strain in (sorted(tuple(strain_list), key=strain_list.get)):
            if strain in cluster_strain_list:
                sql = "update strain_clusters set "
                for i, clust in enumerate(strain_list[strain]):
                    sql = sql + "t" + levels[len(levels) - (i + 1)] + " = " + str(clust) + " ,"
                sql = sql[:-1]
                sql = sql + " where name = '" + strain + "'"
                cur.execute(sql)
                self.snpdb_conn.commit()
            else:
                sql = "insert into strain_clusters ("
                for i, clust in enumerate(strain_list[strain]):
                    sql = sql + "t" + levels[len(levels) - (i + 1)] + ","
                sql = sql + "name ) VALUES ("
                for i, clust in enumerate(strain_list[strain]):
                    sql = sql + str(clust) + ","
                sql = sql + "\'" + strain + "\')"
                cur.execute(sql)
                self.snpdb_conn.commit()
        return strain_list

    def get_clusters(self, cur, co):

        cluster_dict = {}
        cluster_strain_list = {}
        sql = "SELECT name, "

        for cuts in sorted(co, reverse=True, key=int):
            sql = sql + "t" + cuts + " ,"
        sql = sql[:-1]
        sql = sql + " from strain_clusters"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            cluster_strain_list[row[0]] = []
            for i, h in enumerate(sorted(co, key=int)):
                cluster_strain_list[row[0]].append(int(row[i + 1]))
                if i not in cluster_dict:
                    cluster_dict[i] = {}
                    if h not in cluster_dict[i]:
                        cluster_dict[i][row[i + 1]] = []
                        cluster_dict[i][row[i + 1]].append(row[0])
                    else:
                        cluster_dict[i][row[i + 1]].append(row[0])
                elif h not in cluster_dict[i]:
                    cluster_dict[i][row[i + 1]] = []
                    cluster_dict[i][row[i + 1]].append(row[0])
                else:
                    cluster_dict[i][row[i + 1]].append(row[0])

        return cluster_strain_list, cluster_dict

    def merged_clusters(self, cluster_strain_list, strain_list, cur):
        for strain in cluster_strain_list:
            if tuple(cluster_strain_list[strain]) != tuple(strain_list[strain]):
                sql = "insert into cluster_logs (name, old_cluster, new_cluster, date) VALUES (\'%s\',\'%s\' ,\'%s\' ,\'%s\')" % (strain, str(tuple(cluster_strain_list[strain])), str(tuple(strain_list[strain])), str(datetime.now()))
                cur.execute(sql)
                self.snpdb_conn.commit()
                print "!    cluster merge " + strain + " " + str(tuple(cluster_strain_list[strain])) + " to " + str(tuple(strain_list[strain]))

    def update_clusters(self, cur):
        cur.execute('select * from strain_clusters')
        row = cur.fetchall()
        if not row:
            '''
            run the clustering for the first time and add to the db
            '''
            print 'strain_clusters empty'
            print "###  Fetching Matrix:" + str(datetime.now())
            profile_dict = self.get_input(cur)
            cluster_co = self.get_cutoffs(cur)
            clean_clusters = {}
            for cuts in cluster_co:
                links = self.make_links(profile_dict, cuts)
                clusters = self.define_clusters(links)
                clean_clusters[cuts] = self.remove_duplicate_clusters(clusters)
            # self.print_slv_clusters(clean_clusters, cluster_co)
            self.add_clusters_to_table(cur, clean_clusters, cluster_co)
        else:
            '''
            run update_clusters_db
            '''
            profile_dict = self.get_input(cur)
            cluster_co = self.get_cutoffs(cur)
            cluster_strain_list, cluster_dict = self.get_clusters(cur, cluster_co)
            clean_clusters = {}
            for i, cuts in enumerate(sorted(cluster_co, reverse=True, key=int)):
                print "###  Cluster level " + str(cuts) + " :" + str(datetime.now())
                links = self.make_links_update(profile_dict, cuts, cluster_dict[i])
                clusters = self.define_clusters(links)
                clean_clusters[cuts] = self.remove_duplicate_clusters_update(clusters)

            strain_list = self.add_clusters_to_existing_table(clean_clusters, profile_dict, cluster_co, cur, cluster_strain_list)
            self.merged_clusters(cluster_strain_list, strain_list, cur)

class Variant:
    def __init__(self):
        self._pos = int
        self._id = int
        self._product = str
        self._ref_base = str
        self._var_base = str
        self._amino_acid = str
        self._gene = str
        self._locus_tag = str

def vcf_to_db(args, config_dict, vcf):
    logger = logging.getLogger('snapperdb.vcf_to_db')
    logger.info('Initialising SNPdb class')
    snpdb = SNPdb(config_dict)
    logger.info('Parsing config dict')
    snpdb.parse_config_dict(config_dict)
    snpdb._write_conn_string()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)
    if inspect.stack()[0][3] == 'fastq_to_db':
        ## need to add the code to add the single vcf passed to this function to the db.
        pass
    elif inspect.stack()[0][3] == 'vcf_to_db':
        ## there is no existing vcf class here, but there will definitely be a vcf, and there may be a pickle.
        logger.info('Initialising Vcf class')
        vcf = Vcf()
        logger.info('Making SNPdb variables and output files')
        snpdb.define_class_variables_and_make_output_files(args, vcf)
        if os.path.exists(os.path.join(vcf.tmp_dir, vcf.sample_name + '_bad_pos.pick')):
            logger.info('There are already serialised variants and ignored positions for this sample')
            print os.path.join(vcf.tmp_dir, vcf.sample_name + '_bad_pos.pick')
            logger.info('Loading serialised variants and ignored positions')
            bad_pos = pickle.load(open(os.path.join(vcf.tmp_dir, vcf.sample_name + '_bad_pos.pick')))
            good_var = pickle.load(open(os.path.join(vcf.tmp_dir, vcf.sample_name + '_good_var.pick')))
            vcf.good_var = good_var
            vcf.bad_pos = bad_pos
            logger.info('Checking the length of the VCF')
            snpdb.check_len_vcf(vcf)
            logger.info('Uploading to SNPdb')
            snpdb.snpdb_upload(vcf)
        else:
            logger.info('There are no serialised variants, parsing config dict')
            vcf.parse_config_dict(config_dict)
            logger.info('Reading vcf')
            vcf.read_vcf()
            logger.info('Checking length of vcf')
            snpdb.check_len_vcf(vcf)
            logger.info('Serialising variants and ignored positions')
            vcf.pickle_variants_and_ignored_pos(args)
            logger.info('Uploading to SNPdb')
            snpdb.snpdb_upload(vcf)

def make_snpdb(config_dict):
    snpdb = SNPdb(config_dict)
    snpdb._write_conn_string()
    snpdb.make_snpdb()

def read_file(file_name):
    try:
        openfile = open(file_name, 'r')
    except:
        print file_name + " not found ... "
        sys.exit()
    strain_list = []
    for line in openfile:
        strain_list.append(line.strip())
    return strain_list

def read_fasta(ref):
    try:
        openfile = open(ref, 'r')
    except:
        print (ref + " not found ... ")
        sys.exit()
    ref_seq = []
    for line in openfile:
        matchObj = re.search('>', line)
        if matchObj is None:
            for n in line.strip():
                ref_seq.append(n)
    return ref_seq

def read_rec_file(rec_file):
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
    return rec_list

def get_the_snps(args, config_dict):
    print "###  START: " + str(datetime.now())
    print "###  Inititialising SnpDB Class:" + str(datetime.now())
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    strain_list = read_file(args.strain_list)
    snpdb._write_conn_string()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)
    cur = snpdb.snpdb_conn.cursor()
    ref_seq_file = os.path.join(snapperdb.__ref_genome_dir__, snpdb.reference_genome + '.fa')
    ref_seq = read_fasta(ref_seq_file)
    if args.rec_file != 'N':
        print "###  Reading recombination list:" + str(datetime.now())
        rec_list = read_rec_file(args.rec_file)
    else:
        rec_list = []
    snpdb.parse_args_for_get_the_snps(args, cur, strain_list, ref_seq)
    print "###  Printing FASTA: " + str(datetime.now())
    snpdb.print_fasta(args.out, args.alignment_type, args.ref_flag, rec_list)
    if args.mat_flag == 'Y':
        print "###  Printing Matrix: " + str(datetime.now())
        snpdb.print_matrix(args.out)
    if args.var_flag == 'Y':
        print "###  Printing Variants: " + str(datetime.now())
        snpdb.print_vars(args.out, args.alignment_type, rec_list, args.ref_flag)

def update_distance_matrix(config_dict, args):
    print "###  Inititialising SnpDB Class:" + str(datetime.now())
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._write_conn_string()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)
    cur = snpdb.snpdb_conn.cursor()
    print '### Getting strains ' + str(datetime.now())
    strain_list, update_strain = snpdb.get_strains(cur)
    # # get_all_good_ids from snpdb2 takes a snp cutoff as well, here, we don't have a SNP cutoff so we set it arbitrarily high.
    snp_co = '1000000'
    print "###  Populating distance matrix: " + str(datetime.now())
    snpdb.parse_args_for_update_matrix(cur, snp_co, strain_list)
    if args.hpc == 'N':
        print '### Launching serial update_distance_matrix ' + str(datetime.now())
        snpdb.check_matrix(cur, strain_list, update_strain)
        snpdb.update_clusters(cur)
    else:
        try:
            print '### Launching parallele update_distance_matrix ' + str(datetime.now())
            args.hpc = int(args.hpc)
            short_strain_list = set(strain_list) - set(update_strain)
            snpdb.write_qsubs_to_check_matrix(args, strain_list, short_strain_list, update_strain, config_dict['snpdb_name'])
            # # on cluster version this will have to be subject to a qsub hold - no it wont, can just run on headnode
            snpdb.check_matrix(cur, update_strain, update_strain)
        except ValueError as e:
            print '\n#### Error ####'
            print e, '-m has to be an integer'

def qsub_to_check_matrix(config_dict, args):
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._write_conn_string()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)
    cur = snpdb.snpdb_conn.cursor()
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
    snpdb.parse_args_for_update_matrix(cur, snp_co, strain_list)
    snpdb.check_matrix(cur, short_strain_list, update_strain)

    # # need to clean up as otherwise the glob
    os.system('rm -f {0}'.format(args.strain_list))
    direc, name = os.path.split(args.strain_list)
    list_number = name.split('_')[-1]
    shell_script = '{0}/update_matrix_{1}.sh'.format(direc, list_number)
    os.system('rm -f {0}'.format(shell_script))

def update_clusters(config_dict):
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._write_conn_string()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)
    cur = snpdb.snpdb_conn.cursor()
    snpdb.update_clusters(cur)

