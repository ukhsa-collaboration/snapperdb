__author__ = 'flashton'

import psycopg2, psycopg2.extras
import sys
import os
from Bio import SeqIO
import datetime
import inspect
import snapperdb
import errno
from gbru_vcf import Vcf
import pickle
import re
import numpy as np

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

    ## functions below here are for querying the snpdb
    
    def get_background(self,cur,strain_list, args):
        sql = "select distinct("+args.back_flag+") from strain_clusters"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            sql2 = "select sc.name from strain_clusters sc , strain_meta sm where sm.name = sc.name and "+args.back_flag+"="+str(
                row[0])+" "
            if args.meta_flag != 'N':
                temp = args.meta_flag.split(',')
                for meta in temp:
                    temp2 = meta.split(':')
                    sql2 = sql2 + "and "+temp2[0]+"=\'"+str(temp2[1])+"\' "
                sql2 = sql2 + " limit 1"    
            #print sql2
            cur.execute(sql2)
            rows2 = cur.fetchone()
            if rows2:
                strain_list.append(rows2[0])
        return strain_list

    def add_strains_to_sql_co(self, sql, strain_list, co, name):
        for strain in strain_list:
            sql += "name =\'" + strain + "\' or "
        sql = sql[:-4]
        sql += ") and icount("+name+") < "+co
        return sql

    def get_all_good_ids(self,cur,strain_list, snp_co):
        strain_snps = {}
        totlist = []
        ## this could be replaced by where like any (array[strain_list]) - I think
        sql = " select variants_id, name, icount(variants_id) as count from strains_snps where ("
        sql = self.add_strains_to_sql_co(sql, strain_list, snp_co, "variants_id")
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            totlist = set(totlist) | set(row[0])
            strain_snps[row[1]] = row[0]
        return totlist, strain_snps

    def get_variants(self,cur):
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

    def make_consensus(self,ref_seq, ref_flag, cur):
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
                if self.variants[ids].var_base != fasta[strain][self.variants[ids].pos-1]:
                    fasta[strain][self.variants[ids].pos-1] = self.variants[ids].var_base
                    if self.variants[ids].pos:
                        if self.variants[ids].pos in var_look:
                            var_look[self.variants[ids].pos] = var_look[self.variants[ids].pos] + 1
                            if ids not in var_id_list:                            
                                var_id_list.append(ids)
                        else:
                            var_look[self.variants[ids].pos] = 1
                            var_id_list.append(ids)
            for bad_ids in strains_ig:
                fasta[strain][bad_ids-1] = 'N'
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
                        if self.fasta[strain1][var-1] != self.fasta[strain2][var-1] and self.fasta[strain1][var-1] != 'N' and self.fasta[strain2][var-1] !='N':                                 diff_matrix[strain1][strain2]+=1

        return diff_matrix    

    def parse_args_for_get_the_snps(self, args, cur, strain_list, ref_seq):
        ## this populates class variables specific to querying the SNPdb
        if args.back_flag != 'N':
            print "###  Getting background strains:"+str(datetime.datetime.now())
            strain_list = self.get_background(cur, strain_list, args)

        print "###  Getting good positions:"+str(datetime.datetime.now())
        self.goodids, self.strains_snps = self.get_all_good_ids(cur, strain_list, args.snp_co)

        print "Variable positions: "+ str(len(self.goodids))

        if len(self.goodids) == 0:
            print "###  No variable positions found: EXITING "+str(datetime.datetime.now())
            sys.exit()
        print str(len(self.strains_snps)) +" strains used out of "+str(len(strain_list))
        print "###  Getting variants:"+str(datetime.datetime.now())
        self.variants, self.posIDMap = self.get_variants(cur)
        print "###  Making consensus:"+str(datetime.datetime.now())
        self.fasta, self.var_look, self.n_look, self.badlist, self.var_id_list = self.make_consensus(ref_seq, args.ref_flag, cur)
        print "Ignored positions: "+ str(len(self.badlist))
        if args.mat_flag == 'Y':
            print "###  Making Distance Matrix: "+str(datetime.datetime.now())
            self.matrix = self.calc_matrix()

    def print_fasta(self, out, flag, ref_flag, rec_list):
        f = open(out+'.fa','w')
        for strain in self.fasta:
            f.write(">" + strain + "\n")
            if flag == 'W':
                for i, seq in enumerate(self.fasta[strain]):
                    f.write(seq)
            elif flag == 'A':
                for i in sorted(self.var_look):
                    if i not in rec_list:
                        if ref_flag == 'Y':
                            f.write(self.fasta[strain][i-1])
                        elif self.var_look[i] != len(self.strains_snps):
                            f.write(self.fasta[strain][i-1])
            elif flag == 'C':
                for i in sorted(self.var_look):
                    if i not in self.n_look:
                        if i not in rec_list:
                            if ref_flag == 'Y':
                                f.write(self.fasta[strain][i-1])
                            elif self.var_look[i] != len(self.strains_snps):
                                f.write(self.fasta[strain][i-1])
            f.write("\n")

    def print_matrix(self, out):
        f = open(out+'.matrix','w')
        for strain1 in self.matrix:
            for strain2 in self.matrix[strain1]:
                if strain1 != strain2:
                    f.write(strain1+"\t"+strain2+"\t"+str(self.matrix[strain1][strain2])+"\n")
                    
    def print_vars(self, out, flag, rec_list, ref_flag):
        f = open(out+'.variants','w')
        for pos in sorted(self.var_look):
            var_ids = self.posIDMap[pos]
            for var_id in var_ids:
                if var_id in self.var_id_list:
                    if flag == 'W':
                        f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base)  + "\t" + str(self.variants[var_id].gene) + "\n")
                    elif flag == 'A':
                        if pos not in rec_list:
                            if ref_flag == 'Y':
                                f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base)  + "\t" + str(self.variants[var_id].amino_acid)  + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")                    
                            elif self.var_look[pos] != len(self.strains_snps):
                                f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base)  + "\t" + str(self.variants[var_id].amino_acid)  + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                    
                    elif flag == 'C':
                        if pos not in self.n_look:
                            if pos not in rec_list:
                                if ref_flag == 'Y':
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base)  + "\t" + str(self.variants[var_id].amino_acid)  + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")                    
                                elif self.var_look[pos] != len(self.strains_snps):
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base)  + "\t" + str(self.variants[var_id].amino_acid)  + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")

    ## functions below here are from update_distance_matrix

    def get_strains(self, cur):
        strain_list = []
        sql = " select name from strain_stats where ignore is NULL"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            strain_list.append(row[0])
        return strain_list

    def parse_args_for_update_matrix(self, cur, snp_co, strain_list):
        ## this populates class variables specific to querying the SNPdb

        print "###  Getting good positions:"+str(datetime.datetime.now())
        self.goodids, self.strains_snps = self.get_all_good_ids(cur, strain_list, snp_co)
        print "###  Getting variants:"+str(datetime.datetime.now())
        self.variants, self.posIDMap = self.get_variants(cur)

    def get_bad_pos_for_strain_update_matrix(self, cur, strain):
        strain_ig = []
        sql = "select ignored_pos, icount(ignored_pos), name from strains_snps where name =\'" + strain + "\'"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            strain_ig = row[0]
        return strain_ig    

    def check_matrix(self, cur, data_list):
        update_strain = []
        for strain1 in data_list:
            sql = "select * from dist_matrix where strain1 = '%s' or strain2 = '%s' limit 1" % (strain1, strain1)
            cur.execute(sql)
            row = cur.fetchall()
            if not row:
                update_strain.append(strain1)
        seen_strain = []
        for strain1 in update_strain:
            seen_strain.append(strain1)
            print "Populating matrix for: "+ strain1
            strain1_good_var = self.strains_snps[strain1]
            strain1_ig_pos = self.get_bad_pos_for_strain_update_matrix(cur, strain1)
            for strain2 in data_list:
                if strain1 != strain2 and strain2 not in seen_strain:
                    strain2_good_var = self.strains_snps[strain2]
                    strain2_ig_pos = self.get_bad_pos_for_strain_update_matrix(cur, strain2)
                    #getunion of bad_pos
                    all_bad_pos = set(strain1_ig_pos) | set(strain2_ig_pos)

                    #getsymmetric difference of variants
                    all_var = set(strain1_good_var) ^ set(strain2_good_var)
                    diff = 0
                    for var_id in all_var:
                        if self.variants[var_id].pos not in all_bad_pos:
                            diff = diff + 1
                    #add to db
                    sql2 = "insert into dist_matrix (strain1, strain2, snp_dist) VALUES (\'%s\',\'%s\',%s)" % (strain1, strain2, diff)
                    cur.execute(sql2)    
                    self.snpdb_conn.commit()


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
        matchObj = re.search('>',line)
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
            rec_range = range((int(temp[0])-1), (int(temp[1])-1))
            rec_list = set(rec_list) | set(rec_range)
    return rec_list

def get_the_snps(args, config_dict):
    print "###  START: "+str(datetime.datetime.now())
    print "###  Inititialising SnpDB Class:"+str(datetime.datetime.now())
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    strain_list = read_file(args.strain_list)
    snpdb._write_conn_string()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)
    cur = snpdb.snpdb_conn.cursor()
    ref_seq_file = os.path.join(snapperdb.__ref_genome_dir__, snpdb.reference_genome + '.fa')
    ref_seq = read_fasta(ref_seq_file)
    if args.rec_file != 'N':
        print "###  Reading recombination list:"+ str(datetime.datetime.now())
        rec_list = read_rec_file(args.rec_file)
    else:
        rec_list = []
    snpdb.parse_args_for_get_the_snps(args, cur, strain_list, ref_seq)
    print "###  Printing FASTA: "+str(datetime.datetime.now())
    snpdb.print_fasta(args.out, args.alignment_type, args.ref_flag, rec_list)
    if args.mat_flag == 'Y':
        print "###  Printing Matrix: "+str(datetime.datetime.now())
        snpdb.print_matrix(args.out)
    if args.var_flag == 'Y':
        print "###  Printing Variants: "+str(datetime.datetime.now())
        snpdb.print_vars(args.out, args.alignment_type, rec_list, args.ref_flag)

def update_distance_matrix(config_dict):
    print "###  Inititialising SnpDB Class:"+str(datetime.datetime.now())
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._write_conn_string()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)
    cur = snpdb.snpdb_conn.cursor()
    strain_list = snpdb.get_strains(cur)
    ## get_all_good_ids from snpdb2 takes a snp cutoff as well, here, we don't have a SNP cutoff so we set it arbitrarily high.
    snp_co = '1000000'
    print "###  Populating distance matrix:"+str(datetime.datetime.now())
    snpdb.parse_args_for_update_matrix(cur, snp_co, strain_list)
    snpdb.check_matrix(cur, strain_list)