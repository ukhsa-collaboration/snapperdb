__author__ = 'flashton'

from datetime import datetime
import errno
import glob
import os
import sys
from Bio import SeqIO
import psycopg2
import psycopg2.extras
import logging
from variant import Variant
import snapperdb


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
    # # conscious design decision to make cutoffs part of the SNPdb module rather than vcf, as should be consistent within a db
    depth_cutoff = None
    mq_cutoff = None
    ad_cutoff = None
    total_av_depth_co = None

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
            if attr == 'average_depth_cutoff':
                self.average_depth_cutoff = config_dict[attr]

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

    def _connect_to_snpdb(self):
        self.conn_string = 'host=\'{0}\' dbname={1} user=\'{2}\' password=\'{3}\''.format(self.pg_host, self.snpdb_name,
                                                                                          self.pg_uname, self.pg_pword)
        does_snpdb_exist = self._check_if_snpdb_exists()
        if does_snpdb_exist == True:
            self.snpdb_conn = psycopg2.connect(self.conn_string)
        else:
            print 'Cant find snpdb %s' % self.snpdb_name
            pass

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

    def check_duplicate(self, vcf, database):
        dup = False
        dict_cursor = self.snpdb_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        dict_cursor.execute("select distinct(name) FROM %s where name = \'%s\'" % (database, vcf.sample_name))
        for row in dict_cursor:
            dup = True
        dict_cursor.close()
        return dup

    def add_info_to_strain_stats(self, vcf):
        if not self.check_duplicate(vcf, 'strain_stats'):
            time_now = datetime.now()
            time_now = str(time_now)
            insert_statement = 'INSERT INTO strain_stats (name, av_cov, time_of_upload, number_mixed_positions, mixed_positions) ' \
                               'VALUES (%s, %s, %s, %s, %s)'
            cur = self.snpdb_conn.cursor()
            cur.execute(insert_statement, (vcf.sample_name, vcf.depth_average, time_now, vcf.number_mixed_positions,
                                           vcf.mixed_positions))
            self.snpdb_conn.commit()
            cur.close()
        elif self.check_duplicate(vcf, 'strain_stats'):
            sys.stderr.write('%s is already in SNPdb strain_stats %s\n' % (vcf.sample_name, self.reference_genome))

    def add_info_to_strain_stats_mc(self, vcf, mixed_pos_dict):
        if self.check_duplicate(vcf, 'strain_stats') == False:
            mixed_pos_list = []
            for contig in mixed_pos_dict:
                for mixed_pos in mixed_pos_dict[contig]:
                    mixed_pos_list.append('%s position %s' % (contig, mixed_pos))
            time_now = str(datetime.now())
            insert_statement = 'INSERT INTO strain_stats (name, av_cov, time_of_upload, number_mixed_positions, mixed_positions) ' \
                               'VALUES (%s, %s, %s, %s, %s)'
            cur = self.snpdb_conn.cursor()
            cur.execute(insert_statement, (vcf.sample_name, vcf.depth_average, time_now,
                                           vcf.number_mixed_positions, mixed_pos_list))
            self.snpdb_conn.commit()
            cur.close()

    def add_new_variants(self, pos, contig, good_var):
        var_base = good_var[pos]
        # fasta starts at 0, snpdb pos is worked out from a reference that starts at 1 (position in snpdb)
        ref_base = contig[int(pos) - 1]
        cursor = self.snpdb_conn.cursor()
        cursor.execute("insert into variants (pos, var_base, ref_base) VALUES (%s,\'%s\',\'%s\')" % (pos, var_base,
                                                                                                     ref_base))
        self.snpdb_conn.commit()

        cursor.execute('select currval(\'variants_id_seq\')')
        seq_id = cursor.fetchall()
        seq_id = seq_id[0][0]
        return seq_id

    def add_new_variants_mc(self, pos, contig_seq, var_base, contig, cursor):
        ref_base = contig_seq[int(pos) -1]
        cursor.execute("insert into variants (pos, var_base, ref_base, contig) VALUES (%s,\'%s\',\'%s\',\'%s\')" %
                   (pos, var_base, ref_base, contig))
        # print cursor.execute("select currval(\'variants_id_seq\')")

        cursor.execute("select currval(\'variants_id_seq\')")
        res = cursor.fetchall()
        seq_id = str(res[0][0])
        return int(seq_id)

    def add_new_ig_pos_mc(self, cursor, contig, pos):
        cursor.execute('insert into ignored_pos (pos, contig) VALUES (%s, \'%s\')' % (pos, contig))
        cursor.execute("select currval(\'ignored_pos_id_seq\')")
        res = cursor.fetchall()
        seq_id = res[0][0]
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

    def add_to_snpdb_mc(self, bad_pos_dict, var_dict, mixed_pos_dict, sample_name):
        ref_genome_path = os.path.join(self.ref_genome_dir, self.reference_genome + '.fa')
        ref_fasta_dict = {}
        if os.path.exists(ref_genome_path):
            with open(ref_genome_path, 'r') as fi:
                ref_fasta = SeqIO.parse(fi, 'fasta')
                for contig in ref_fasta:
                    ref_fasta_dict[contig.id] = contig.seq
        else:
            sys.stderr.write('Could not find {0}, check your file extension (needs to be .fa)\n'.format
                             (ref_genome_path))
            sys.exit()
        dict_cursor = self.snpdb_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        cursor = self.snpdb_conn.cursor()
        var_db_list = []
        for contig in var_dict:
            existing_variants_dict = {}
            query = 'SELECT * FROM variants where contig = %s'
            queryObj = dict_cursor.execute(query, (contig,))
            if queryObj != None:
                res = queryObj.dictresult()
                for row in res:
                    if row['pos'] in existing_variants_dict:
                        existing_variants_dict[row['pos']][row['var_base']] = row['id']
                    else:
                        existing_variants_dict[row['pos']] = {}
                        existing_variants_dict[row['pos']][row['var_base']] = row['id']

            ## on the adding new variants bit
            contig_seq = ref_fasta_dict[contig]
            for pos in var_dict[contig]:
                if int(pos) not in existing_variants_dict:
                    seq_id = self.add_new_variants_mc(pos, contig_seq, var_dict[contig][pos], contig, cursor)
                    var_db_list.append(seq_id)
                elif var_dict[contig][pos] not in existing_variants_dict[int(pos)]:
                    seq_id = self.add_new_variants_mc(pos, contig_seq, var_dict[contig][pos], contig, cursor)
                    var_db_list.append(seq_id)
                else:
                    var_db_list.append(existing_variants_dict[int(pos)][var_dict[contig][pos]])

        ig_db_list = []
        for contig in bad_pos_dict:
            contig_seq = ref_fasta_dict[contig]
            ig_dic = {}
            query = 'SELECT * FROM ignored_pos where contig = %s'
            queryObj = dict_cursor.execute(query, (contig,))
            if queryObj != None:
                res = queryObj.dictresult()
                for row in res:
                    if row['pos'] in ig_dic:
                        ig_dic[row['pos']][row['var_base']] = row['id']
                    else:
                        ig_dic[row['pos']] = {}
                        ig_dic[row['pos']][row['var_base']] = row['id']

            for pos in bad_pos_dict[contig]:
                if pos not in ig_dic:
                    seq_id = self.add_new_ig_pos_mc(cursor, contig, pos)
                    ig_db_list.append(seq_id)
                else:
                    ig_db_list.append(ig_dic[int(pos)])

        insert_statement = 'insert into strains_snps (name, variants_id, ignored_pos) values (%s, %s, %s)'
        cursor.execute(insert_statement, (sample_name, var_db_list, ig_db_list))
        self.snpdb_conn.commit()

    def snpdb_upload(self, vcf):
        if not self.check_duplicate(vcf, 'strains_snps'):
            self.add_info_to_strain_stats(vcf)
            print 'depth is', vcf.depth_average, self.average_depth_cutoff
            if vcf.depth_average >= int(self.average_depth_cutoff):
                self.add_to_snpdb(vcf)
            else:
                update_statement = 'UPDATE strain_stats SET ignore = \'i - average depth below cutoff\' where name = \'%s\' ' \
                                   % vcf.sample_name
                cur = self.snpdb_conn.cursor()
                cur.execute(update_statement)
                self.snpdb_conn.commit()
                cur.close()
                sys.stderr.write('average depth below cutoff, not added to SNPdb')
        elif self.check_duplicate(vcf, 'strains_snps'):
            sys.stderr.write('%s is already in SNPdb strains_snps %s\n' % (vcf.sample_name, self.reference_genome))

    def snpdb_upload_multi_contig(self, vcf, bad_pos_dict, var_dict, mixed_pos_dict):
        if not self.check_duplicate(vcf, 'strains_snps'):
            self.add_info_to_strain_stats_mc(vcf, mixed_pos_dict)
            if vcf.depth_average >= int(self.average_depth_cutoff):
                self.add_to_snpdb_mc(bad_pos_dict, var_dict, mixed_pos_dict, vcf.sample_name)

    # # functions below here are for querying the snpdb

    def get_background(self, strain_list, args):
        cur = self.snpdb_conn.cursor()
        sql = 'SELECT name FROM strain_clusters WHERE id IN (SELECT MIN(id) FROM strain_clusters GROUP BY %s)' % args.back_flag
        # sql = "select distinct(%s) from strain_clusters"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            strain_list.append(row[0])

        ## not currently handling meta - perhaps should add?

        return strain_list

    def add_strains_to_sql_co(self, sql, strain_list, co, name):
        for strain in strain_list:
            sql += "name =\'" + strain + "\' or "
        sql = sql[:-4]
        sql += ") and icount(" + name + ") < " + co
        return sql

    def get_all_good_ids(self, strain_list, snp_co):
        cur = self.snpdb_conn.cursor()
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

    def get_variants(self):
        cur = self.snpdb_conn.cursor()
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

    def get_variants_mc(self):
        cur = self.snpdb_conn.cursor()
        variant_container = {}
        pos_2_id_list = {}
        sql = "select pos , id, ref_base, var_base, contig from variants"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            variant = Variant()
            variant.pos = row[0]
            variant.id = row[1]
            variant.ref_base = row[2]
            variant.var_base = row[3]
            variant.contig = row[4]
            variant_container[row[1]] = variant
            if variant.id in self.goodids:
                if row[4] in pos_2_id_list:
                    if row[0] in pos_2_id_list[row[4]]:
                        pos_2_id_list[row[4]][row[0]].append(row[1])
                    else:
                        pos_2_id_list[row[4]][row[0]] = []
                        pos_2_id_list[row[4]][row[0]].append(row[1])
                else:
                    pos_2_id_list[row[4]] = {}
                    pos_2_id_list[row[4]][row[0]] = []
                    pos_2_id_list[row[4]][row[0]].append(row[1])
        return variant_container, pos_2_id_list

    def get_bad_pos_for_strain_get_the_vars(self, strain, totlist):
        cur = self.snpdb_conn.cursor()
        strain_ig = {}
        sql = "select ignored_pos, icount(ignored_pos), name as count from strains_snps where name =\'" + strain + "\'"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            totlist = set(totlist) | set(row[0])
            strain_ig = row[0]
        return totlist, strain_ig

    def make_consensus(self, ref_seq, ref_flag):
        fasta = {}
        var_id_list = []
        var_look = {}
        n_look = {}
        badlist = []
        for strain in self.strains_snps:
            badlist, strains_ig = self.get_bad_pos_for_strain_get_the_vars(strain, badlist)
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

    def make_consensus_mc(self, ref_seq, args, reference_genome_name):
        fasta = {}
        var_look = {}
        n_look = {}
        for strain in self.strains_snps:
            if strain not in fasta:
                fasta[strain] = dict(ref_seq)
            for ids in sorted(self.strains_snps[strain]):
                if self.variants[ids].var_base != fasta[strain][self.variants[ids].contig][self.variants[ids].pos-1] \
                        and ids not in self.strains_snps[reference_genome_name]:
                    fasta[strain][self.variants[ids].contig][self.variants[ids].pos-1] = self.variants[ids].var_base
                    if self.variants[ids].pos-1:
                        if self.variants[ids].contig not in var_look:
                            var_look[self.variants[ids].contig] = {}
                            var_look[self.variants[ids].contig][self.variants[ids].pos-1] = 'V'
                        else:
                            var_look[self.variants[ids].contig][self.variants[ids].pos-1] = 'V'

            for bad_ids in self.igpos[strain]:
                fasta[strain][self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos-1] = 'N'
                if self.IgPos_container[bad_ids].contig not in n_look:
                    n_look[self.IgPos_container[bad_ids].contig] = {}
                    n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos-1] = 'N'
                else:
                    n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos-1] = 'N'
            for bad_ids in self.igpos[reference_genome_name]:
                fasta[strain][self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos-1] = 'N'
                if self.IgPos_container[bad_ids].contig not in n_look:
                    n_look[self.IgPos_container[bad_ids].contig] = {}
                    n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos-1] = 'N'
                else:
                    n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos-1] = 'N'

        # fasta['ref'] = deepcopy(ref_seq)
        return fasta, var_look, n_look

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

    def parse_args_for_get_the_snps(self, args, strain_list, ref_seq):
        logger = logging.getLogger('snapperdb.SNPdb.parse_args_for_get_the_snps')

        # # this populates class variables specific to querying the SNPdb
        if args.back_flag != 'N':
            logger.info('Getting background strains')
            strain_list = self.get_background(strain_list, args)
        logger.info('Getting good positions')
        self.goodids, self.strains_snps = self.get_all_good_ids(strain_list, args.snp_co)
        logger.info('Variable positions: ' + str(len(self.goodids)))
        if len(self.goodids) == 0:
            logger.error('No variable positions found: EXITING')
            sys.exit()
        logger.info(str(len(self.strains_snps)) + ' strains used out of ' + str(len(strain_list)))
        logger.info('Getting variants')
        self.variants, self.posIDMap = self.get_variants()
        logger.info('Making consensus')
        self.fasta, self.var_look, self.n_look, self.badlist, self.var_id_list = self.make_consensus(ref_seq,
                                                                                                     args.ref_flag)
        logger.info('Ignored positions' + str(len(self.badlist)))
        if args.mat_flag == 'Y':
            logger.info('Making Distance Matrix')
            self.matrix = self.calc_matrix()

    def parse_args_for_get_the_snps_mc(self, args, strain_list, ref_seq, reference_genome_name):
        self.goodids, self.strains_snps = self.get_all_good_ids(strain_list, args.snp_co)

        if len(self.goodids) == 0:
            print 'No variable positions found: EXITING'
            # logger.error('No variable positions found: EXITING')
            sys.exit()
        print (str(len(self.strains_snps)) + ' strains used out of ' + str(len(strain_list)))
        self.variants, self.posIDMap = self.get_variants_mc()
        self.fasta, self.var_look, self.n_look, self.badlist, self.var_id_list = \
            self.make_consensus_mc(ref_seq, args, reference_genome_name)




    def print_fasta(self, out, flag, rec_list, ref_flag):
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

    def get_strains(self):
        cur = self.snpdb_conn.cursor()
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

    def parse_args_for_update_matrix(self, snp_co, strain_list):
        # # this populates class variables specific to querying the SNPdb
        cur = self.snpdb_conn.cursor()
        print "###  Getting good positions:" + str(datetime.now())
        self.goodids, self.strains_snps = self.get_all_good_ids(strain_list, snp_co)
        print "###  Getting variants:" + str(datetime.now())
        self.variants, self.posIDMap = self.get_variants()

    def get_bad_pos_for_strain_update_matrix(self, strain):
        cur = self.snpdb_conn.cursor()
        strain_ig = []
        sql = "select ignored_pos, icount(ignored_pos), name from strains_snps where name =\'" + strain + "\'"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            strain_ig = row[0]
        return strain_ig

    def check_matrix(self, data_list, update_strain):
        cur = self.snpdb_conn.cursor()
        seen_strain = []
        for strain1 in update_strain:
            seen_strain.append(strain1)
            print "Populating matrix for: " + strain1
            strain1_good_var = self.strains_snps[strain1]
            strain1_ig_pos = self.get_bad_pos_for_strain_update_matrix(strain1)
            for strain2 in data_list:
                if strain1 != strain2 and strain2 not in seen_strain:
                    strain2_good_var = self.strains_snps[strain2]
                    strain2_ig_pos = self.get_bad_pos_for_strain_update_matrix(strain2)
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
        logs_dir = '/phengs/hpc_projects/routine_salmonella_snapperdb/logs_chron_jobs'
        scripts_dir = '/phengs/hpc_projects/routine_salmonella_snapperdb/scripts_chron_job/non_fastq_to_vcf'
        this_dir = os.path.dirname(os.path.realpath(__file__))
        snapperdb_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

        # os.system('rm -rf {0}/update_list*'.format(this_dir))
        # os.system('rm -rf {0}/update_matrix_*'.format(this_dir))


        for i, each in enumerate(self.chunks(update_strain, args.hpc)):
            with open('{0}/{1}.{3}.update_list_{2}'.format(scripts_dir, args.now, i, self.snpdb_name), 'w') as fo:
                for x in each:
                    fo.write(x + '\n')
        with open('{0}/{1}.{2}.short_strain_list'.format(scripts_dir, args.now, self.snpdb_name), 'w') as fo:
            for x in short_strain_list:
                fo.write(x + '\n')

        with open('{0}/{1}.{2}.strain_list'.format(scripts_dir, args.now, self.snpdb_name), 'w') as fo:
            for x in strain_list:
                fo.write(x + '\n')

        res = sorted(glob.glob('{0}/{1}.{2}.update_list*'.format(scripts_dir, args.now, self.snpdb_name)))

        for i, update_list in enumerate(res):
            command = ('#! /bin/bash\n'
                       '#$ -o {0}/{7}.check_matrix.stdout\n'
                       '#$ -e {0}/{7}.check_matrix.stderr\n'
                       '#$ -m e\n'
                       '#$ -wd {1}\n'
                       '#$ -N up_mat_{2}_{3}\n\n'
                       '. /etc/profile.d/modules.sh\n'
                       'module load snapperdb/0.1\n'
                       'SnapperDB_main.py'
                       ' qsub_to_check_matrix -c {4}'
                       ' -l {5}/{7}.{8}.strain_list'
                       ' -s {5}/{7}.{8}.short_strain_list'
                       ' -u {6}\n'.format(logs_dir, logs_dir, snpdb, i, args.config_file, scripts_dir, update_list, args.now, self.snpdb_name))

            with open('{0}/{2}.{3}.update_mat_{1}.sh'.format(scripts_dir, i, args.now, self.snpdb_name), 'w') as fo:
                fo.write(command)
            os.system('qsub {0}/{2}.{3}.update_mat_{1}.sh'.format(scripts_dir, i, args.now, self.snpdb_name))


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

    def get_input(self):
        cur = self.snpdb_conn.cursor()
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

    def get_cutoffs(self):
        cur = self.snpdb_conn.cursor()
        co = []
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

    def add_clusters_to_table(self, clusters, levels):
        cur = self.snpdb_conn.cursor()
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

    def add_clusters_to_existing_table(self, clusters, profile_dict, levels, cluster_strain_list):
        cur = self.snpdb_conn.cursor()
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

    def get_clusters(self, co):
        cur = self.snpdb_conn.cursor()
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

    def merged_clusters(self, cluster_strain_list, strain_list):
        cur = self.snpdb_conn.cursor()
        for strain in cluster_strain_list:
            if tuple(cluster_strain_list[strain]) != tuple(strain_list[strain]):
                sql = "insert into cluster_logs (name, old_cluster, new_cluster, date) VALUES (\'%s\',\'%s\' ,\'%s\' ,\'%s\')" % (strain, str(tuple(cluster_strain_list[strain])), str(tuple(strain_list[strain])), str(datetime.now()))
                cur.execute(sql)
                self.snpdb_conn.commit()
                print "!    cluster merge " + strain + " " + str(tuple(cluster_strain_list[strain])) + " to " + str(tuple(strain_list[strain]))

    def update_clusters(self):
        cur = self.snpdb_conn.cursor()
        cur.execute('select * from strain_clusters')
        row = cur.fetchall()
        if not row:
            '''
            run the clustering for the first time and add to the db
            '''
            print 'strain_clusters empty'
            print "###  Fetching Matrix:" + str(datetime.now())
            profile_dict = self.get_input()
            cluster_co = self.get_cutoffs()
            clean_clusters = {}
            for cuts in cluster_co:
                links = self.make_links(profile_dict, cuts)
                clusters = self.define_clusters(links)
                clean_clusters[cuts] = self.remove_duplicate_clusters(clusters)
            # self.print_slv_clusters(clean_clusters, cluster_co)
            self.add_clusters_to_table(clean_clusters, cluster_co)
        else:
            '''
            run update_clusters_db
            '''
            profile_dict = self.get_input()
            cluster_co = self.get_cutoffs()
            cluster_strain_list, cluster_dict = self.get_clusters(cluster_co)
            clean_clusters = {}
            for i, cuts in enumerate(sorted(cluster_co, reverse=True, key=int)):
                print "###  Cluster level " + str(cuts) + " :" + str(datetime.now())
                links = self.make_links_update(profile_dict, cuts, cluster_dict[i])
                clusters = self.define_clusters(links)
                clean_clusters[cuts] = self.remove_duplicate_clusters_update(clusters)

            strain_list = self.add_clusters_to_existing_table(clean_clusters, profile_dict, cluster_co, cluster_strain_list)
            self.merged_clusters(cluster_strain_list, strain_list, )

    ## functions below here are for getting variants of interest

    def x(self):
        '''
        To do
        '''
        pass
