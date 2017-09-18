__author__ = 'gidis'

from datetime import datetime
import errno
import glob
import os
import math
import json
import sys
import re
from Bio import SeqIO
import psycopg2
import psycopg2.extras
import logging
from variant import Variant
import snapperdb
from copy import deepcopy

class Igpos:
    def __init__(self):
        self._pos = int
        self._id = int
        self._contig = str
# -------------------------------------------------------------------------------------------------

class SNPdb:
    
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
# -------------------------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------------------------

    def mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise
# -------------------------------------------------------------------------------------------------

    def make_tmp_dir(self, args):
        self.tmp_dir = os.path.join(os.path.dirname(args.vcf[0]), 'snpdb')
        self.mkdir_p(self.tmp_dir)

# -------------------------------------------------------------------------------------------------

    def define_class_variables_and_make_output_files(self, args, vcf):
        # # need to handle either vcf or fastqs
        try:
            vcf.sample_name = os.path.basename(args.vcf[0]).split(os.extsep)[0]
        except AttributeError:
            vcf.sample_name = os.path.basename(args.fastqs[0]).split(os.extsep)[0]

        vcf.ref_genome_path = os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fa')
        vcf.make_tmp_dir(args)

        try:
            vcf.vcf_filehandle = args.vcf[0]
        except:
            vcf.vcf_filehandle = os.path.join(vcf.tmp_dir, '{0}.filtered.vcf'.format(vcf.sample_name))

# -------------------------------------------------------------------------------------------------

    def _connect_to_snpdb(self):
        self.conn_string = 'host=\'{0}\' dbname={1} user=\'{2}\' password=\'{3}\''.format(self.pg_host, self.snpdb_name,
                                                                                          self.pg_uname, self.pg_pword)
        does_snpdb_exist = self.check_if_snpdb_exists()
        if does_snpdb_exist == True:
            self.snpdb_conn = psycopg2.connect(self.conn_string)
        else:
            print '### Cant connect to SnapperDB %s' % self.snpdb_name
            #print '### Please check database exists and Postgres server is running'
            #sys.exit()

# -------------------------------------------------------------------------------------------------

    def check_if_snpdb_exists(self):
        try:
            psycopg2.connect(self.conn_string)
            return True
        except psycopg2.OperationalError:
            return False

# -------------------------------------------------------------------------------------------------
    def make_snpdb(self):
        does_snpdb_exist = self.check_if_snpdb_exists()
        if does_snpdb_exist == True:
            sys.stderr.write(self.snpdb_name + ' already exists\n')
        else:
            sys.stdout.write('The SNPdb {0} does not exist - running sql to create database\n'.format(self.snpdb_name))
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
            sql_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'template_snapperdb_denovo_refs_sql')
            cur.execute(open(sql_script, 'r').read())
            conn.commit()
            conn.close()
        return does_snpdb_exist

# -------------------------------------------------------------------------------------------------

    def add_cluster(self):
        cur = self.snpdb_conn.cursor()
        cur.execute("insert into strain_clusters (name, t250, t100, t50, t25, t10, t5, t0) VALUES (\'%s\',1,1,1,1,1,1,1)" % self.reference_genome)
        self.snpdb_conn.commit()

 # -------------------------------------------------------------------------------------------------

    def check_duplicate(self, vcf, database):
        dup = False
        dict_cursor = self.snpdb_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        dict_cursor.execute("select distinct(name) FROM %s where name = \'%s\'" % (database, vcf.sample_name))
        for row in dict_cursor:
            dup = True
        dict_cursor.close()
        return dup

# -------------------------------------------------------------------------------------------------

    def add_info_to_strain_stats(self, vcf):
        if self.check_duplicate(vcf, 'strain_stats') == False:
            time_now = datetime.now()
            time_now = str(time_now)
            insert_statement = 'INSERT INTO strain_stats (name, av_cov, time_of_upload, number_mixed_positions) ' \
                               'VALUES (%s, %s, %s, %s)'
            cur = self.snpdb_conn.cursor()
            cur.execute(insert_statement, (vcf.sample_name, vcf.depth_average, time_now, vcf.number_mixed_positions))
            self.snpdb_conn.commit()
            cur.close()

# -------------------------------------------------------------------------------------------------


    def add_new_variants(self, pos, ref_base, var_base, contig, cursor):
        cursor.execute("insert into variants (pos, var_base, ref_base, contig) VALUES (%s,\'%s\',\'%s\',\'%s\')" %
                   (pos, var_base, ref_base, contig))
        cursor.execute("select currval(\'variants_id_seq\')")
        res = cursor.fetchall()
        seq_id = str(res[0][0])
        return int(seq_id)

# -------------------------------------------------------------------------------------------------

    def add_new_ig_pos(self, cursor, contig, pos):
        cursor.execute('insert into ignored_pos (pos, contig) VALUES (%s, \'%s\')' % (pos, contig))
        cursor.execute("select currval(\'ignored_pos_id_seq\')")
        res = cursor.fetchall()
        seq_id = res[0][0]
        return seq_id

# -------------------------------------------------------------------------------------------------

    def snpdb_annotate_vars(self, vcf):
        #get genbank


        ref_gbk_path = os.path.join(self.ref_genome_dir, self.reference_genome + '.gbk')
        if os.path.exists(ref_gbk_path):
            #get variant objects that arnet annotated
            self.variants = self.get_variants_annotate()
            self.get_gbk(ref_gbk_path)
            self.add_annotate_vars_to_db()

# -------------------------------------------------------------------------------------------------

    def add_annotate_vars_to_db(self):
        for variant_id in self.variants:
            cur = self.snpdb_conn.cursor()
            sql = "update variants set locus_tag = %s, amino_acid = %s, product = %s, gene = %s where id = %s"
            cur.execute(sql, (self.variants[variant_id].locus_tag, self.variants[variant_id].amino_acid, self.variants[variant_id].product, self.variants[variant_id].gene, variant_id))
            self.snpdb_conn.commit()
            cur.close()

  # -------------------------------------------------------------------------------------------------


    def get_gbk(self, gb_file):
        try:
            gb_records = list(SeqIO.parse(gb_file, "genbank"))
            for variant_id in self.variants:
                flag = 0
                pos = self.variants[variant_id].pos
                for gb_record in gb_records:
                    if gb_record.id in self.variants[variant_id].contig:
                        for feature in gb_record.features:
                            if feature.type == 'CDS':
                                if pos >= feature.location.start and pos < feature.location.end:
                                        flag = 1
                                        ref_sequence = gb_record.seq[feature.location.start:feature.location.end]
                                        var_sequence = gb_record.seq[feature.location.start:pos-1] + self.variants[variant_id].var_base + gb_record.seq[pos:feature.location.end]
                                        ref_prot = ""
                                        var_prot = ""
                                        if feature.strand == 1:
                                            ref_prot = ref_sequence.translate(table=11)
                                            var_prot = var_sequence.translate(table=11)
                                        elif feature.strand == -1:
                                            ref_prot = ref_sequence.reverse_complement().translate(table=11)
                                            var_prot = var_sequence.reverse_complement().translate(table=11)
                                        if str(ref_prot) ==  str(var_prot):
                                            self.variants[variant_id].amino_acid = 'SYNONYMOUS'
                                        else:
                                            for i, base in enumerate(ref_prot):
                                                if ref_prot[i] != var_prot[i]:
                                                    self.variants[variant_id].amino_acid = str(ref_prot[i]) + str(i+1) + str(var_prot[i])
                                        if "product" in feature.qualifiers:
                                            self.variants[variant_id].product = str(feature.qualifiers['product'][0])
                                        else:
                                            self.variants[variant_id].product = ''
                                        if "gene" in feature.qualifiers:
                                            self.variants[variant_id].gene = str(feature.qualifiers['gene'][0])
                                        else:
                                            self.variants[variant_id].gene = ''
                                        if "locus_tag" in feature.qualifiers:
                                            self.variants[variant_id].locus_tag = str(feature.qualifiers['locus_tag'][0])
                                        else:
                                            self.variants[variant_id].locus_tag = ''

                if flag == 0:
                    self.variants[variant_id].amino_acid = 'NON CODING'
                    self.variants[variant_id].product = ''
                    self.variants[variant_id].gene = ''
                    self.variants[variant_id].locus_tag = ''
        except IOError:
            print 'Cannot find {0}'.format(gb_file)
# -------------------------------------------------------------------------------------------------


    def query_snpdb(self,vcf):
        print "Under development"


# -------------------------------------------------------------------------------------------------

    def add_to_snpdb(self, vcf):
        #get reference genome
        ref_genome_path = os.path.join(self.ref_genome_dir, self.reference_genome + '.fa')
        ref_fasta_dict = {}
        #open reference genome
        if os.path.exists(ref_genome_path):
            with open(ref_genome_path, 'r') as fi:
                ref_fasta = SeqIO.parse(fi, 'fasta')
                # get each of the contigs
                for contig in ref_fasta:
                    ref_fasta_dict[contig.id] = contig.seq
        else:
            #CHANGE this needs to be logged and and try and except
            sys.stderr.write('### Could not find {0}, check your file extension (needs to be .fa)\n'.format
                             (ref_genome_path))
            sys.exit()

        #get a cursor
        dict_cursor = self.snpdb_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        cursor = self.snpdb_conn.cursor()
        #initialise a list for vars
        var_db_list = []
        ig_db_list = []


        #loop through the different contigs in vcf
        for parsed_vcf in vcf.parsed_vcf_container:
            #get any known variants in that contig
            existing_variants_dict = {}
            query = 'SELECT * FROM variants where contig = %s'
            dict_cursor.execute(query, (parsed_vcf.ref,))
            res = dict_cursor.fetchall()
            if len(res) != 0:
                 for row in res:
                    try:
                        existing_variants_dict[row['pos']][row['var_base']] = row['id']
                    except:
                        existing_variants_dict[row['pos']] = {}
                        existing_variants_dict[row['pos']][row['var_base']] = row['id']

            #go through the variants in this contig
            for pos in parsed_vcf.good_var:
                #if we havent seen a variant at this poistion before
                if int(pos) not in existing_variants_dict:
                    seq_id = self.add_new_variants(pos, parsed_vcf.ref_base[pos], parsed_vcf.var[pos], parsed_vcf.ref, cursor)
                    var_db_list.append(seq_id)
                #if we have seen a variant at this position but not rhis variant
                elif parsed_vcf.var[pos] not in existing_variants_dict[int(pos)]:
                    seq_id = self.add_new_variants(pos, parsed_vcf.ref_base[pos], parsed_vcf.var[pos], parsed_vcf.ref, cursor)
                    var_db_list.append(seq_id)
                else:
                    var_db_list.append(existing_variants_dict[int(pos)][parsed_vcf.var[pos]])


            #get any ignored_positions know variants in that contig
            ig_dic = {}
            query = 'SELECT * FROM ignored_pos where contig = %s'
            dict_cursor.execute(query, (parsed_vcf.ref,))
            res = dict_cursor.fetchall()
            if len(res) != 0:
                for row in res:
                    try:
                        ig_dic[row['pos']] = row['id']
                    except:
                        ig_dic[row['pos']] = {}
                        ig_dic[row['pos']] = row['id']

            #go through the ignored_pos in this contig
            for pos in parsed_vcf.bad_pos:
                if int(pos) not in ig_dic:
                    seq_id = self.add_new_ig_pos(cursor, parsed_vcf.ref, pos)
                    ig_db_list.append(seq_id)
                else:
                    ig_db_list.append(ig_dic[int(pos)])

        #add into strains_snp
        insert_statement = 'insert into strains_snps (name, variants_id, ignored_pos) values (%s, %s, %s)'
        cursor.execute(insert_statement, (vcf.sample_name, var_db_list, ig_db_list))
        self.snpdb_conn.commit()
# -------------------------------------------------------------------------------------------------

    def is_number(self,s):
        try:
            float(s)
            return True
        except ValueError:
            return False

# -------------------------------------------------------------------------------------------------

    def snpdb_upload(self, vcf,args):
        #lets check if they are forcing
        
        if not self.check_duplicate(vcf, 'strains_snps'):
            print 'Calulated depth is %s - cuttoff is %s' % (vcf.depth_average, self.average_depth_cutoff)
            if args.force == 'Y' or (self.is_number(vcf.depth_average) == True and  float(vcf.depth_average) >= float(self.average_depth_cutoff)):
                #add strains
                if (args.force == 'Y' and self.is_number(vcf.depth_average) == False):
                    vcf.depth_average = None

                self.add_info_to_strain_stats(vcf)
                #CHANGE - this needs to be logged
                self.add_to_snpdb(vcf)
            else:
                update_statement = 'UPDATE strain_stats SET ignore = \'i - average depth below cutoff\' where name = \'%s\' ' \
                                   % vcf.sample_name
                cur = self.snpdb_conn.cursor()
                cur.execute(update_statement)
                self.snpdb_conn.commit()
                cur.close()
                #CHANGE this needs to be logged
                sys.stderr.write('Average depth below cutoff, not added to SNPdb\n')
        elif self.check_duplicate(vcf, 'strains_snps'):
            #CHANGE this needs to be logged
            sys.stderr.write('%s is already in SNPdb strains_snps %s\n' % (vcf.sample_name, self.reference_genome))

# -------------------------------------------------------------------------------------------------


    def snpdb_query(self, vcf,args):
        
        #need to return something like a strains_snps dict
        self.query_snpdb(vcf)


    # # functions below here are for querying the snpdb

# -------------------------------------------------------------------------------------------------


    def get_background(self, strain_list, args):
        cur = self.snpdb_conn.cursor()
        sql = 'SELECT name FROM strain_clusters WHERE id IN (SELECT MIN(id) FROM strain_clusters GROUP BY %s)' % args.back_flag
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            strain_list.append(row[0])

        return strain_list
# -------------------------------------------------------------------------------------------------

    def add_strains_to_sql_co(self, sql, strain_list, co, name):
        for strain in strain_list:
            sql += "name =\'" + strain + "\' or "
        sql = sql[:-4]
        sql += ") and icount(" + name + ") < " + co
        return sql
# -------------------------------------------------------------------------------------------------

    def get_all_good_ids(self, strain_list, snp_co):
        #set up cursor
        cur = self.snpdb_conn.cursor()
        #dict of strain_name : list sof variants_ids
        strain_snps = {}
        # list of all variants_ui
        totlist = []
        # build sql
        sql = " select variants_id, name, icount(variants_id) as count from strains_snps where ("
        sql = self.add_strains_to_sql_co(sql, strain_list, snp_co, "variants_id")
        cur.execute(sql)
        rows = cur.fetchall()
        #populate data structures
        for row in rows:
            totlist = set(totlist) | set(row[0])
            strain_snps[row[1]] = row[0]
        return totlist, strain_snps

# -------------------------------------------------------------------------------------------------

    def get_variants_mc(self):
        cur = self.snpdb_conn.cursor()
        #set up dict for Variant objects where key is variant id
        variant_container = {}
        #set up dict of contig : pos : id
        pos_2_id_list = {}
        sql = "select pos , id, ref_base, var_base, contig, amino_acid, gene, product from variants"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            variant = Variant()
            variant.pos = row[0]
            variant.id = row[1]
            variant.ref_base = row[2]
            variant.var_base = row[3]
            variant.contig = row[4]
            variant.amino_acid = row[5]
            variant.gene = row[6]
            variant.product = row[7]
            variant_container[row[1]] = variant
            # if this variant is present in the list of strains we are interested in add to pos_2_id_list
            
            if variant.id in self.goodids:
                try:
                     pos_2_id_list[row[4]][row[0]].append(row[1])
                except:
                    try:
                            pos_2_id_list[row[4]][row[0]] = []
                            pos_2_id_list[row[4]][row[0]].append(row[1])
                    except:
                            pos_2_id_list[row[4]] = {}
                            pos_2_id_list[row[4]][row[0]] = []
                            pos_2_id_list[row[4]][row[0]].append(row[1])                                       
        return variant_container, pos_2_id_list
# -------------------------------------------------------------------------------------------------

    def get_variants_annotate(self):
        cur = self.snpdb_conn.cursor()
        #set up dict for Variant objects where key is variant id
        variant_container = {}
        #set up dict of contig : pos : id
        pos_2_id_list = {}
        sql = "select pos , id, ref_base, var_base, contig from variants where amino_acid is NULL order by contig, pos"
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
            # if this variant is present in the list of strains we are interested in add to pos_2_id_list
        return variant_container

# -------------------------------------------------------------------------------------------------

    def get_bad_pos_mc(self):
        cur = self.snpdb_conn.cursor()
        #list of all bad positons
        totlist = []
        #dict of bad positions per strains
        ig_pos = {}
        query = 'select ignored_pos, icount(ignored_pos), name from strains_snps where name in %s'
        strain_names = tuple(self.strains_snps.keys())
        cur.execute(query, (strain_names,))
        res = cur.fetchall()
        for row in res:
            totlist = set(totlist) | set(row[0])
            ig_pos[row[2]] = set(row[0])
        return totlist, ig_pos

# -------------------------------------------------------------------------------------------------

    def get_igs_mc(self):
        cur = self.snpdb_conn.cursor()
        #create ignored pos container where key is id
        igPos_container = {}
        pos_2_id_list = {}

        sql = "select pos, id, contig from ignored_pos"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            ig_pos = Igpos()
            ig_pos.pos = row[0]
            ig_pos.id = row[1]
            ig_pos.contig = row[2]
            igPos_container[row[1]] = ig_pos

            try:
                pos_2_id_list[row[2]][row[0]] = row[1]
            except:
                pos_2_id_list[row[2]] = {}
                pos_2_id_list[row[2]][row[0]] = row[1]                            

        return igPos_container, pos_2_id_list

# -------------------------------------------------------------------------------------------------

    def make_consensus_mc(self, ref_seq, args, reference_genome_name):
        #create a  dictionary for the alignment
        fasta = {}
        #create a dict to capture variants
        var_look = {}
        #create a dict to capture ignored pos
        n_look = {}
        #list of variants to output
        var_id_list = set()

        #for strain
        for strain in self.strains_snps:
            if strain != reference_genome_name:
                # if we haven't seen this strain create a deepcopy to mutate
                if strain not in fasta:
                    fasta[strain] = deepcopy(ref_seq)

                # go the contig ids
                for ids in sorted(self.strains_snps[strain]):
                    #if this base is not the same as the reference and not a variant when the reference is mapped against iteself
                    if self.variants[ids].var_base != fasta[strain][self.variants[ids].contig][self.variants[ids].pos-1] \
                            and ids not in self.strains_snps[reference_genome_name]:
                        #change the reference position to the variant base
                        fasta[strain][self.variants[ids].contig][self.variants[ids].pos-1] = self.variants[ids].var_base
                        #capture the number of variants per position
                        if self.variants[ids].pos:
                            try:
                                var_look[self.variants[ids].contig][self.variants[ids].pos] = var_look[self.variants[ids].contig][self.variants[ids].pos] + 1
                            except:
                                try:
                                    var_look[self.variants[ids].contig][self.variants[ids].pos] = 1
                                except:
                                    var_look[self.variants[ids].contig] = {}
                                    var_look[self.variants[ids].contig][self.variants[ids].pos] = 1
                # go through ignored_pos
                for bad_ids in self.igpos[strain]:
                    #set the reference position to a N
                    fasta[strain][self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos-1] = 'N'
                    try:
                        n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] = n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] + 1
                    except:
                        try:
                            n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] = 1
                        except:
                            n_look[self.IgPos_container[bad_ids].contig] = {}
                            n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] = 1

                # add the ignored pos from the referecne genome
                for bad_ids in self.igpos[reference_genome_name]:
                    fasta[strain][self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos-1] = 'N'
                    try:
                        n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] = n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] + 1
                    except:
                        try:
                            n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] = 1
                        except:
                            n_look[self.IgPos_container[bad_ids].contig] = {}
                            n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] = 1                           

        # create a deepcopy to get the reference back
        if args.ref_flag == 'Y':
            fasta[reference_genome_name] = deepcopy(ref_seq)
             # add the ignored pos back in from the referecne genome
            for bad_ids in self.igpos[reference_genome_name]:
                fasta[reference_genome_name][self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos-1] = 'N'
                try:
                    n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] = n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] + 1
                except:
                    try:
                        n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] = 1
                    except:
                        n_look[self.IgPos_container[bad_ids].contig] = {}
                        n_look[self.IgPos_container[bad_ids].contig][self.IgPos_container[bad_ids].pos] = 1     
        else:
              del self.strains_snps[reference_genome_name]
        return fasta, var_look, n_look,var_id_list
# -------------------------------------------------------------------------------------------------


    def calc_matrix_mc(self):
        diff_matrix = {}
        for strain1 in self.fasta:
            diff_matrix[strain1] = {}
            for strain2 in self.fasta:
                if strain2 in diff_matrix and strain1 not in diff_matrix[strain2]:
                    diff_matrix[strain1][strain2] = 0
                    for contig in self.var_look:
                        for var in self.var_look[contig]:
                            if self.fasta[strain1][contig][var-1] != self.fasta[strain2][contig][var-1] and self.fasta[strain1][contig][var-1] != 'N' and self.fasta[strain2][contig][var-1] !='N':
                                diff_matrix[strain1][strain2]+=1
        return diff_matrix
# -------------------------------------------------------------------------------------------------


    def parse_args_for_get_the_snps_mc(self, args, strain_list, ref_seq, reference_genome_name):
        logger = logging.getLogger('snapperdb.SNPdb.parse_args_for_get_the_snps')

        #get some background strains if requested
        if args.back_flag != 'N':
            logger.info('Getting background strains')
            strain_list = self.get_background(strain_list, args)

        #get all varitants and variants for each strains
        logger.info('Getting good positions')
        print '### Revcovering Variants...'
        self.goodids, self.strains_snps = self.get_all_good_ids(strain_list, args.snp_co)
        print '### Variable positions: ' + str(len(self.goodids))
        logger.info('Variable positions: ' + str(len(self.goodids)))

        #If there are no variants returned we can exit
        if len(self.goodids) == 0:
            print '### No variable positions found: EXITING'
            sys.exit()

        print '### '+str(len(self.strains_snps)) + ' strains used out of ' + str(len(strain_list))
        logger.info(str(len(self.strains_snps)) + ' strains used out of ' + str(len(strain_list)))

        print '### Recovering ignored positions...'
        logger.info('Getting ignored positions')
        self.badlist, self.igpos = self.get_bad_pos_mc()
        print '### Ignored positions: ' + str(len(self.badlist))
        logger.info('Ignored positions: ' + str(len(self.badlist)))

        #get actual variant objects
        print '### Creating Variant Objects'
        self.variants, self.posIDMap = self.get_variants_mc()
        # get ignored position objects
        print '### Creating Ignored Positions Objects'

        self.all_bad_pos = []
        self.IgPos_container, self.igposIDMap = self.get_igs_mc()
        print '### Creating Consensus... '

        # make consensus sequence based on the above
        self.fasta, self.var_look, self.n_look, self.var_id_list = self.make_consensus_mc(ref_seq, args, reference_genome_name)
        if args.mat_flag == 'Y':
            self.matrix = self.calc_matrix_mc()
# -------------------------------------------------------------------------------------------------


    def print_fasta_mc(self,args, rec_dict):
        #open fasta
        f = open(args.out + '.fa', 'w')
        for strain in self.fasta:
            #write header
            f.write(">" + strain + "\n")
            #for each contig
            for contig in sorted(self.fasta[strain]):
                #if flag is for a whole genome alignment add whole contig
                if args.alignment_type == 'W':
                    for i, seq in enumerate(self.fasta[strain][contig]):
                        f.write(seq)
                #if we want to include Ns in the alignment
                elif args.alignment_type == 'A':
                    if contig in self.var_look:
                        for i in sorted(self.var_look[contig]):
                            #get number of N's for this position
                            ncount = 0
                            if contig in self.n_look:
                                if i in self.n_look[contig]:
                                    ncount = self.n_look[contig][i]
                            #check recombination
                            rec_flag = False
                            if contig in rec_dict:
                                if i in rec_dict[contig]:
                                    rec_flag = True
                            if rec_flag == False:
                                if args.ref_flag == 'Y':
                                    #print as long as all no ref posotions not an N
                                    if i in self.n_look[contig]:
                                        if self.n_look[contig][i] < (len(self.strains_snps)-1):
                                            f.write(self.fasta[strain][contig][i-1])
                                    else:
                                        f.write(self.fasta[strain][contig][i-1])
                                #print as long as not all var or N
                                elif self.var_look[contig][i]+ncount < len(self.strains_snps):
                                    f.write(self.fasta[strain][contig][i-1])
                #if we want to include a proportion of N's
                elif re.match("A:(\d+)" ,args.alignment_type) is not None:
                    m = re.match("A:(\d+)" ,args.alignment_type)
                    co = m.group(1)
                    if contig in self.var_look:
                        for i in sorted(self.var_look[contig]):
                            #get number of N's for this position
                            ncount = 0
                            if contig in self.n_look:
                                if i in self.n_look[contig]:
                                    ncount = self.n_look[contig][i]
                            #check recombination
                            rec_flag = False
                            if contig in rec_dict:
                                if i in rec_dict[contig]:
                                    rec_flag = True
                            if rec_flag == False:
                                if args.ref_flag == 'Y':
                                    #print as long as all no ref posotions not an N
                                    if i in self.n_look[contig]:
                                        if (float(float(self.n_look[contig][i]) / float(len(self.strains_snps)-1)*100) < float(100-float(co))):
                                            f.write(self.fasta[strain][contig][i-1])
                                    else:
                                            f.write(self.fasta[strain][contig][i-1])
                                elif self.var_look[contig][i]+ncount < len(self.strains_snps):
                                    if i in self.n_look[contig]:
                                        #dont print N if this posisition in all N
                                        if (float(float(self.n_look[contig][i]) / float(len(self.strains_snps))*100) < float(100-float(co))):
                                            f.write(self.fasta[strain][contig][i-1])
                                    else:
                                            f.write(self.fasta[strain][contig][i-1])
                else:
                    if contig in self.var_look:
                        for i in sorted(self.var_look[contig]):
                            if contig in rec_dict:
                                if i not in rec_dict[contig]:
                                    if contig not in self.n_look:
                                        f.write(self.fasta[strain][contig][i-1])
                                    elif i not in self.n_look[contig]:
                                        f.write(self.fasta[strain][contig][i-1])
                            else:
                                if contig not in self.n_look:
                                        f.write(self.fasta[strain][contig][i-1])
                                elif i not in self.n_look[contig]:
                                    f.write(self.fasta[strain][contig][i-1])
            f.write("\n")


# -------------------------------------------------------------------------------------------------

    def print_matrix(self, out):
        f = open(out + '.matrix', 'w')
        for strain1 in self.matrix:
            for strain2 in self.matrix[strain1]:
                if strain1 != strain2:
                    f.write(strain1 + "\t" + strain2 + "\t" + str(self.matrix[strain1][strain2]) + "\n")

# -------------------------------------------------------------------------------------------------


    def print_vars_mc(self,args, rec_dict):
        f = open(args.out + '.variants', 'w')
        for contig in sorted(self.var_look):
            for pos in sorted(self.var_look[contig]):
                var_ids = self.posIDMap[contig][pos]
                for var_id in var_ids:
                    if args.alignment_type == 'W':
                        f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].gene) + "\n")
                    elif args.alignment_type == 'A':
                       #get number of N's for this position
                        ncount = 0
                        if contig in self.n_look:
                            if pos in self.n_look[contig]:
                                ncount = self.n_look[contig][pos]

                        #check recombination
                        rec_flag = False
                        if contig in rec_dict:
                            if pos in rec_dict[contig]:
                                rec_flag = True
                        if rec_flag == False:
                            if args.ref_flag == 'Y':
                                if pos in self.n_look[contig]:
                                    #dont print N if this posisition in all N
                                    if self.n_look[contig][pos] < (len(self.strains_snps)-1):
                                         f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                                else:
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                            #This is to remove positions that are all variants from the reference genoe
                            elif self.var_look[contig][pos]+ncount < len(self.strains_snps):
                                if pos in self.n_look[contig]:
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                                else:
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                    elif re.match("A:(\d+)" ,args.alignment_type) is not None:
                        m = re.match("A:(\d+)" ,args.alignment_type)
                        co = m.group(1)
                        #get number of N's for this position
                        ncount = 0
                        if contig in self.n_look:
                            if pos in self.n_look[contig]:
                                ncount = self.n_look[contig][pos]
                        #check recombination
                        rec_flag = False
                        if contig in rec_dict:
                            if pos in rec_dict[contig]:
                                rec_flag = True
                        if rec_flag == False:
                            if args.ref_flag == 'Y':
                                if pos in self.n_look[contig]:
                                    #dont print N if this posisition in all N
                                    if (float(float(self.n_look[contig][pos]) / float(len(self.strains_snps)-1)*100) < float(100-float(co))):
                                         f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                                else:
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                            #This is to remove positions that are all variants from the reference genoe
                            elif self.var_look[contig][pos]+ncount < len(self.strains_snps):
                                #dont print N if this posisition in all N
                                if pos in self.n_look[contig]:
                                    if (float(float(self.n_look[contig][pos]) / float(len(self.strains_snps))*100) < float(100-float(co))):
                                        f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                                else:
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                    else:
                        if pos not in self.n_look[contig]:
                            if contig in rec_dict:
                                if pos not in rec_dict[contig]:
                                    if args.ref_flag == 'Y':
                                        f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                                    elif self.var_look[contig][pos] != len(self.strains_snps):
                                        f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                            else:
                                if args.ref_flag == 'Y':
                                        f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")
                                elif self.var_look[contig][pos] != len(self.strains_snps):
                                    f.write(str(var_id) + "\t" + str(self.variants[var_id].contig) + "\t" + str(self.variants[var_id].pos) + "\t" + str(self.variants[var_id].var_base) + "\t" + str(self.variants[var_id].amino_acid) + "\t" + str(self.variants[var_id].gene) + "\t" + str(self.variants[var_id].product) + "\n")


    # # functions below here are from update_distance_matrix

    def get_strains(self):
        cur = self.snpdb_conn.cursor()
        strain_list = []
        sql = "select name from strain_stats where ignore is NULL"
        cur.execute(sql)
        rows = cur.fetchall()
        strain_list = [row[0] for row in rows]

        sql = "select distinct(strain1) from dist_matrix"
        cur.execute(sql)
        rows = cur.fetchall()
        row_strain1 = [row[0] for row in rows]

        sql = "select distinct(strain2) from dist_matrix"
        cur.execute(sql)
        rows = cur.fetchall()
        row_strain2 = [row[0] for row in rows]

        all_strains = set(row_strain1) | set(row_strain2)

        update_strain = set(strain_list) -  set(all_strains)

        return strain_list, update_strain

    def parse_args_for_update_matrix(self, snp_co, strain_list):


        logger = logging.getLogger('snapperdb.SNPdb.parse_args_for_update_matrix')
        logger.info('Getting good positions')
        self.goodids, self.strains_snps = self.get_all_good_ids(strain_list, snp_co)
        logger.info('Variable positions: ' + str(len(self.goodids)))
        #get actual variant objects
        self.variants, self.posIDMap = self.get_variants_mc()
        #get ignored positions object
        self.IgPos_container, self.igposIDMap = self.get_igs_mc()


    def get_bad_pos_for_strain_update_matrix(self, strain):
        cur = self.snpdb_conn.cursor()
        strain_ig = []
        sql = "select ignored_pos, icount(ignored_pos), name from strains_snps where name =\'" + strain + "\'"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            strain_ig = row[0]
        return strain_ig


    def write_qsubs_to_check_matrix(self, args, idx, one_strain, present_strains, snpdb):
        '''
        how to call this particular function as a qsub - need access to check matrix on command line.
        1. write lists
        2. for each in lists, qsub
        '''

        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        logs_dir = os.getcwd()
        scripts_dir = os.getcwd()

        with open('{0}/{1}.{2}.added_strain.txt'.format(scripts_dir, idx, self.snpdb_name), 'w') as fo:
            for x in one_strain:
                fo.write(x + '\n')

        with open('{0}/{1}.{2}.present_strains.txt'.format(scripts_dir, idx, self.snpdb_name), 'w') as fo:
            for x in present_strains:
                fo.write(x + '\n')

        command = ('#! /bin/bash\n'
                       '#$ -o {0}/{6}.check_matrix.stdout\n'
                       '#$ -e {0}/{6}.check_matrix.stderr\n'
                       '#$ -m e\n'
                       '#$ -wd {0}\n'
                       '#$ -N up_mat_{1}_{2}\n\n'
                       '. /etc/profile.d/modules.sh\n'
                       'module purge\n'
                       'module load gastro/snapperdb/0.2.2\n'
                       'module load uge/8.1.5\n'
                       'SnapperDB_main.py'
                       ' qsub_to_check_matrix -c {3}'
                       ' -l {4}/{2}.{6}.added_strain.txt'
                       ' -s {4}/{2}.{6}.present_strains.txt\n'
                       .format(logs_dir, snpdb, idx, args.config_file, scripts_dir, args.now, self.snpdb_name))

        with open('{0}/{2}.{3}.update_mat_{1}.sh'.format(scripts_dir, idx, timestamp, self.snpdb_name), 'w') as fo:
            fo.write(command)
        os.system('qsub {0}/{2}.{3}.update_mat_{1}.sh'.format(scripts_dir, idx, timestamp, self.snpdb_name))





# -------------------------------------------------------------------------------------------------

    def check_matrix(self, data_list, update_strain):
        cur = self.snpdb_conn.cursor()
        seen_strain = set()

        lookup = range(6500000)
        strain_ig_pos_dict = {}
        for strn in set(data_list):
            strain_ig_pos_dict[strn] = [lookup[x] if x < 6500000 else x for x in self.get_bad_pos_for_strain_update_matrix(strn)]

        newrows = []
        #add the reference genomes bad positions
        ref_ig_pos = set()
        ref_ig_pos = self.get_bad_pos_for_strain_update_matrix(self.reference_genome)

        for strain1 in update_strain:
            seen_strain.add(strain1)
            print "Populating matrix for: " + strain1
            # get ids of all variants in strain1
            strain1_good_var = self.strains_snps[strain1]
            # get ids of all bad positions in strain1
            # strain1_ig_pos = self.get_bad_pos_for_strain_update_matrix(strain1)
            strain1_ig_pos = set(strain_ig_pos_dict[strain1])

            # don't loop over the ones already seen
            data_set = set(data_list).difference(seen_strain)

            for strain2 in data_set:
                # get ids of all variants in strain2
                strain2_good_var = self.strains_snps[strain2]
                # get ids of all bad positions in strain2
                # strain2_ig_pos = self.get_bad_pos_for_strain_update_matrix(strain2)
                strain2_ig_pos = strain_ig_pos_dict[strain2]
                # get union of bad position ids
                all_bad_ids = strain1_ig_pos | set(strain2_ig_pos) | set(ref_ig_pos)
                # get symmetric difference of variant ids
                all_var = set(strain1_good_var) ^ set(strain2_good_var)
                # get sets of (position, contig) tuples
                all_var_pos = set([(self.variants[var_id].pos, self.variants[var_id].contig) for var_id in all_var])
                all_bad_pos = set([(self.IgPos_container[bad_id].pos, self.IgPos_container[bad_id].contig) for bad_id in all_bad_ids])
                # the difference is the number of variants that are not at a bad position
                diff = len(all_var_pos.difference(all_bad_pos))
                newrows.append((strain1, strain2, diff))
        # add to db
        sql2 = "insert into dist_matrix (strain1, strain2, snp_dist) VALUES (%s, %s, %s)"
        cur.executemany(sql2, tuple(newrows))
        self.snpdb_conn.commit()

# -------------------------------------------------------------------------------------------------


    def get_input(self):
        cur = self.snpdb_conn.cursor()
        dist_mat = {}
        sql = "select strain1, strain2, snp_dist from dist_matrix"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            try:
                dist_mat[row[0]][row[1]] = row[2]
            except:
                dist_mat[row[0]] = {}
                dist_mat[row[0]][row[1]] = row[2]
            try:
                dist_mat[row[1]][row[0]] = row[2]
            except:
                dist_mat[row[1]] = {}
                dist_mat[row[1]][row[0]] = row[2]
        return dist_mat

# -------------------------------------------------------------------------------------------------


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

# -------------------------------------------------------------------------------------------------


    def make_links(self, dist_mat, co, cluster_dict):
        clusters = {}
        seen_strain = []
        #for each set of clusters
        for i, cluster in enumerate(sorted(cluster_dict, key=int)):
            #take a copy of the clusters
            clusters[cluster] = cluster_dict[cluster]
            # for each strain in the defined cluster
            for strain1 in cluster_dict[cluster]:
                #add it to the seen list
                seen_strain.append(strain1)
                #compare all the strains to this and add to new cluster if within threshold and not already in there
                for strain2 in dist_mat[strain1]:
                    if int(dist_mat[strain1][strain2]) <= int(co) and strain2 not in cluster_dict[cluster]:
                        clusters[cluster].append(strain2)
                        seen_strain.append(strain2)

        #count what number cluster we are up to
        cluster_count = int(cluster)

        #go through the distance matrix
        for i, strain1 in enumerate(dist_mat):
                #if the strain has not been seen before we need to add it to a cluster
                if strain1 not in seen_strain:
                    cluster_count = cluster_count + 1
                    clusters[cluster_count] = []
                    clusters[cluster_count].append(strain1)
                    #compare to all other strains and see if we can join any other strains to this cluster
                    for strain2 in dist_mat[strain1]:
                        if int(dist_mat[strain1][strain2]) <= int(co) and strain2 not in seen_strain:
                            clusters[cluster_count].append(strain2)

        return clusters

# -------------------------------------------------------------------------------------------------


    def define_clusters(self, links):
        clusters = {}
        for each in links:
            clusters[each] = set(links[each])
            for each2 in links:
                    #check to see if the clusters have anysample in common
                    int_set = set(links[each]) & set(links[each2])
                    if len(int_set) >= 1:
                        #if they have take the union between them
                        clusters[each] = set(clusters[each]) | set(links[each2])
                        # loop through again to see if we can merge with any other clusters
                        for made_clusters in clusters:
                            #check to see if the clusters have anysample in common
                            int_set = set(clusters[each]) & set(clusters[made_clusters])
                            if len(int_set) >= 1:
                                #if they have take the union between them
                                clusters[each] = set(clusters[made_clusters]) | set(clusters[each])
                                #update made clusters
                                clusters[made_clusters] = clusters[each]
        return clusters

# -------------------------------------------------------------------------------------------------


    def remove_duplicate_clusters(self, clusters):
        clean_clusters = {}
        #get clusters with largest members
        for cluster1 in clusters:
            if tuple(sorted(clusters[cluster1])) not in clean_clusters:
                clean_clusters[tuple(sorted(clusters[cluster1]))] = cluster1
            elif int(cluster1) < int(clean_clusters[tuple(sorted(clusters[cluster1]))]):
                clean_clusters[tuple(sorted(clusters[cluster1]))] = cluster1
        return clean_clusters

# -------------------------------------------------------------------------------------------------


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

# -------------------------------------------------------------------------------------------------


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

# -------------------------------------------------------------------------------------------------


    def add_clusters_to_existing_table(self, clusters, dist_mat, levels, cluster_strain_list):
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


# -------------------------------------------------------------------------------------------------

    def get_clusters(self):
        cur = self.snpdb_conn.cursor()
        co = []
        cluster_dict = {}
        cluster_strain_list = {}

        #get columns names
        sql = "SELECT * FROM information_schema.columns WHERE table_schema = 'public' AND table_name = 'strain_clusters'"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            if row[3][0] == "t":
                co.append(row[3][1:])

        co = sorted(co,key=int)
        sql = "SELECT name, "
        for cuts in sorted(co, reverse=True,key=int):
            sql = sql + "t" + cuts + " ,"
        sql = sql[:-1]
        sql = sql+ " from strain_clusters"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            cluster_strain_list[row[0]] = []
            for i, h in enumerate(sorted(co, key=int)):
                cluster_strain_list[row[0]].append(int(row[i+1]))
                try:
                    cluster_dict[i][row[i+1]].append(row[0])
                except:
                    try:
                        cluster_dict[i][row[i+1]] = []
                        cluster_dict[i][row[i+1]].append(row[0])
                    except:
                        cluster_dict[i] = {}
                        cluster_dict[i][row[i+1]] = []
                        cluster_dict[i][row[i+1]].append(row[0])
        return co, cluster_strain_list, cluster_dict

# -------------------------------------------------------------------------------------------------


    def merged_clusters(self, cluster_strain_list, strain_list):
        cur = self.snpdb_conn.cursor()
        for strain in cluster_strain_list:
            if tuple(cluster_strain_list[strain]) != tuple(strain_list[strain]):
                sql = "insert into cluster_logs (name, old_cluster, new_cluster, date) VALUES (\'%s\',\'%s\' ,\'%s\' ,\'%s\')" % (strain, str(tuple(cluster_strain_list[strain])), str(tuple(strain_list[strain])), str(datetime.now()))
                cur.execute(sql)
                self.snpdb_conn.commit()
                print "!    cluster merge " + strain + " " + str(tuple(cluster_strain_list[strain])) + " to " + str(tuple(strain_list[strain]))

# -------------------------------------------------------------------------------------------------

    def get_outliers(self):
        cur = self.snpdb_conn.cursor()
        outliers = []
        sql = "select name from strain_stats where zscore_check = 'Y'"
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            outliers.append(row[0])
        return outliers

# -------------------------------------------------------------------------------------------------

    def check_clusters(self,clusters, dist_mat,levels, cluster_strain_list,outliers):
        strain_list = {}
        cluster_list = {}
        bad_list = []
        #define cluster list
        for co in (sorted(clusters, key=clusters.get)):
            cluster_list[co] = {}
            for i, cluster in enumerate(clusters[co]):
                cluster_list[co][clusters[co][cluster]] = cluster
                for strain in list(cluster):
                    try:
                        strain_list[strain].append(int(clusters[co][cluster]))
                    except KeyError:
                        strain_list[strain] = [int(clusters[co][cluster])]

        #for new strains look at zscore for that cluster
        seen_clusters = {}
        for strain in (sorted(tuple(strain_list), key=strain_list.get)):
            if strain not in cluster_strain_list:
                #print new results
                print strain, strain_list[strain]
                for i, level in enumerate(sorted(levels, reverse=True, key=int)):

                    if level in seen_clusters and strain_list[strain][i] in seen_clusters[level]:
                        print "### Already Investigated ",
                        print level,strain_list[strain][i]
                    else:
                        print "### Investigating ",
                        print level,strain_list[strain][i]

                        try:
                            seen_clusters[level].append(strain_list[strain][i])
                        except KeyError:
                            seen_clusters[level] = [strain_list[strain][i]]

                        check_list = set(cluster_list[level][strain_list[strain][i]])

                        total_dist = {}
                        total = 0
                        for strain1 in check_list:

                            strain2_list = set(dist_mat[strain1]).intersection(check_list)
                            total_dist[strain1] = sum([dist_mat[strain1][strain2] for strain2 in strain2_list])
                            total += total_dist[strain1]
                            total_dist[strain1] = float(total_dist[strain1]) / float(len(check_list))

                        av = float(total) / (float(len(check_list))*float(len(check_list)))
                        print "### Average Distance ",
                        print av
                        if av > 0:
                            sd_top = 0
                            for strain1 in total_dist:
                                sd_top = sd_top + ((av-total_dist[strain1])*(av-total_dist[strain1]))
                            sd = math.sqrt(sd_top / (len(check_list)-1))
                            print "### Standard Deviation ",
                            print sd
                            if sd > 0:
                                for strain1 in total_dist:
                                    z_score = (total_dist[strain1]-av)/ sd
                                    if z_score <= -1.75:
                                        if strain1 not in outliers:
                                            bad_list.append(strain1)
                                            print "### Outlier Z-Score ",
                                            print strain1, total_dist[strain1], z_score
        return bad_list

# -------------------------------------------------------------------------------------------------


    def update_clusters(self):
        #get distance matrix
        dist_mat = self.get_input()
        #get clusters
        co, cluster_strain_list, cluster_dict = self.get_clusters()
        #print cluster_dict
        clean_clusters = {}
        for i, cuts in enumerate(sorted(co,reverse=True,key=int)):
            print "### Cluster level "+str(cuts)+" :"+ str(datetime.time(datetime.now()))
            #print "making links"
            links = self.make_links(dist_mat, cuts, cluster_dict[i])
            #print "defining_clusters"
            clusters = self.define_clusters(links)
            #print "removing duplicates"
            clean_clusters[cuts] = self.remove_duplicate_clusters(clusters)

        print "### Getting previously checked outliers:"+ str(datetime.time(datetime.now()))
        outliers = self.get_outliers()

        print "### Checking Clusters:"+ str(datetime.time(datetime.now()))
        bad_list = self.check_clusters(clean_clusters, dist_mat, co,cluster_strain_list,outliers)

        if not bad_list:
            strain_list = self.add_clusters_to_existing_table(clean_clusters, dist_mat, co, cluster_strain_list)
            self.merged_clusters(cluster_strain_list, strain_list)
        else:
            print "### Not Updating Clusters:"+ str(datetime.time(datetime.now()))


# -------------------------------------------------------------------------------------------------

    def sql_single_extrac(self,data, data_type, table):
        cur = self.snpdb_conn.cursor()
        sql = "SELECT "+ data_type + " FROM " + table + " WHERE name = %s"
        cur.execute(sql, (data,))
        value = (cur.fetchone())[0]
        return value        

# -------------------------------------------------------------------------------------------------


    def parse_args_for_export(self, args, strain_list, ref_seq):
        
        logger = logging.getLogger('snapperdb.SNPdb.parse_args_for_export')
        snp_co = '1000000'
   
        self.goodids, self.strains_snps = self.get_all_good_ids(strain_list,snp_co)
        self.badlist, self.igpos = self.get_bad_pos_mc()
        self.variants, self.posIDMap = self.get_variants_mc()
        self.IgPos_container, self.igposIDMap = self.get_igs_mc()


        for strain in strain_list:
            self.output_dict = {}
            print "###  Exporting JSON for "+strain+": " + str(datetime.time(datetime.now()))

            #get strain stats average coverage
            self.output_dict['strain_stats'] = self.sql_single_extrac(strain, "av_cov", "strain_stats")
            #get variants
            self.convert_seq_for_export(ref_seq,strain)
            self.output_dict['sample'] = strain
            self.output_dict['config_file'] = args.config_file

            self.write_json(strain)


  # -------------------------------------------------------------------------------------------------


    def convert_seq_for_export(self,ref_seq,strain):          

        base_tuple = ('A', 'C', 'G', 'T')
        for s_contig in ref_seq:
            pos_dict = {}
            for base in base_tuple:
                pos_list = []
                for var_id in self.strains_snps[strain]:
                    if base in self.variants[var_id].var_base and s_contig in self.variants[var_id].contig:
                        pos_list.append(str(self.variants[var_id].pos)+"."+self.variants[var_id].ref_base)
                    else:
                        continue
                    if not pos_list:
                        pos_dict[base] = [0]
                    else:
                        pos_dict[base] = pos_list
            pos_list = []

            if s_contig in self.igposIDMap:
                for pos in self.igposIDMap[s_contig]:
                    if self.igposIDMap[s_contig][pos] in self.igpos[strain]:
                        pos_list.append(pos)
            if not pos_list:
                pos_dict['N'] = [0]
            else:
                pos_dict['N'] = pos_list
            # adding gaps as place holder for now
            pos_dict['-'] = [0]
            self.output_dict[s_contig] = pos_dict


  # -------------------------------------------------------------------------------------------------

    def write_json(self,strain):
        #translate dictionnary to Json
        output_file = strain + "_json.txt"
        outfile = open(output_file, 'w')
        outfile.write(json.dumps(self.output_dict))
        outfile.close()
        bash_command = "tar -czvf " + output_file + ".tar.gz " + output_file
        os.system(bash_command)

