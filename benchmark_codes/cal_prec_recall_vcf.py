import vcfpy
import logging
import sys
import numpy as np
import sqlite3
import pandas as pd
import os
import uuid
import subprocess
import itertools as it
from pathlib import Path


logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(name)s:%(exc_info)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)
logger.setLevel(logging.INFO)


def create_sqlite_db(dbpath=":memory:", force=False):
    if os.path.exists(dbpath):
        if force:
            os.remove(dbpath)
            Path(dbpath).touch()
            
    conn = None
    try:
        conn = sqlite3.connect(dbpath)
        logger.info("Establish connection with database at {} with sqlite3 version {}".format(dbpath,sqlite3.sqlite_version))
    except Exception as e:
        logger.error(e)
        raise e
    finally:
        return conn
    

def db_recording_benchmark(dbpath=":memory:", table_map={"":""}, upsert_sql=[""], output_tab_path={"":""}):
    conn = create_sqlite_db(dbpath)
    if conn is not None:
        '''
        First make sure the table and relevant fields exists
        '''
        cursor = conn.cursor()
        for tab_name, create_query in table_map.items():
            try:
                cursor.execute(create_query)
            except Exception as e:
                logger.error(e)
                raise e
        conn.commit()
        '''
        Then insert the value specified by sqlite query
        '''
        for sql in upsert_sql:
            try:
                cursor.execute(sql)
            except Exception as e:
                logger.error(e)
                logger.info(sql)
                raise e
        conn.commit()
        cursor.close()
        '''
        Output table to tab file
        '''
        for tab_name, path in output_tab_path.items():
            '''
            Here needs to make the file replacement operation atomic
            '''
            tmp_path = path + "." + str(uuid.uuid4())
            pd.read_sql_query("SELECT * FROM {}".format(tab_name), conn).replace([None], np.nan).to_csv(tmp_path, sep="\t", index=False)
            result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_path, path, path), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            logger.info(result.stdout)
            
        conn.close()
    else:
        logger.error("Failed to establish connection with sqlite3 db: {}".format(dbpath))
    
    
    
def create_sql_tab(tab_name, db_conn=None, pk = "primary_key", keys=[""], **kwargs):
    tab_create_sql = " CREATE TABLE IF NOT EXISTS {} ( ".format(tab_name) + \
                     ", ".join([ "{} {}".format(k,v) for k,v in kwargs.items() ]) + \
                     ", CONSTRAINT {} PRIMARY KEY(".format(pk) + \
                     ",".join(keys) + "));"
                     
    logger.info("Composed this SQL query clause: \n{}".format(tab_create_sql))
    if db_conn is not None:
        cursor = db_conn.cursor()
        try:
            cursor.execute(tab_create_sql)
        except Exception as e:
            logger.error(e)
            raise e
        db_conn.commit()
        cursor.close()
        return db_conn
    else:
        return tab_create_sql
    
    
def compose_mod_sql_query(tabname="", keys=[], **kwargs):
    '''
    kwargs should be pairs that k is the name of var and v is the value of var
    '''
    mod_sql_query = "INSERT INTO {} ( ".format(tabname) + \
                    ", ".join([ k for k,v in kwargs.items() ]) + \
                    ") VALUES (" + \
                    ", ".join([ "\"{}\"".format(v) if type(v) == str and v != "NULL" else "{}".format(v) for k,v in kwargs.items() ]) + \
                    ") ON CONFLICT ({}) DO UPDATE SET ".format(", ".join(keys)) + \
                    ", ".join([ "{}=excluded.{}".format(k, k) for k,v in kwargs.items() ]) + \
                    ";"
    logger.info("Composed this SQL query clause: \n{}".format(mod_sql_query))
    return mod_sql_query
    

def compare_per_pair(called_vcf, golden_vcf, meta_file, rough = False, exonic = False, gatk = False):
    samp_id = os.path.basename(called_vcf).split(".")[0]
    creader = vcfpy.Reader.from_path(called_vcf)
    
    # Prepare the called recs
    raw_crecs = [ r for r in creader ]
    if exonic: 
        raw_crecs = [ r for r in raw_crecs if "exonic" in r.INFO["Func.refGene"] or "splic" in r.INFO["Func.refGene"] or "UTR" in r.INFO["Func.refGene"] ]
    if gatk:
        raw_crecs = [ r for r in raw_crecs if "GATK" in r.FILTER ]
    crecs = list(dict.fromkeys([ vcf_rec(r, rough=rough) for r in raw_crecs if "LIKELY_INTRINSIC" not in r.FILTER and r.calls[0].gt_type in [1,2] ]))
    crecs_num = len(set([ vcf_rec(r, onlysite=True) for r in raw_crecs if "LIKELY_INTRINSIC" not in r.FILTER and r.calls[0].gt_type in [1,2] ]))
    intrinsic_crecs = list(dict.fromkeys([ vcf_rec(r, rough=rough) for r in raw_crecs if "LIKELY_INTRINSIC" in r.FILTER and r.calls[0].gt_type in [1,2] ]))
    intrinsic_crecs_num = len(set([ vcf_rec(r, onlysite=True) for r in raw_crecs if "LIKELY_INTRINSIC" in r.FILTER and r.calls[0].gt_type in [1,2] ]))
    
    # Prepare the gold recs
    greader = vcfpy.Reader.from_path(golden_vcf)
    raw_grecs = [ r for r in greader ]
    if exonic:
        raw_grecs = [ r for r in raw_grecs if "exonic" in r.INFO["Func.refGene"] or "splic" in r.INFO["Func.refGene"] ]
    grecs = [ vcf_rec(r, rough=rough) for r in raw_grecs ]
    grecs_num = len(grecs)
    
    try:
        recall = len([ r for r in grecs if r in crecs ])/len(grecs)
    except ZeroDivisionError:
        recall = "NULL"
        
    try:
        prec = len([r for r in crecs if r in grecs])/len(crecs)
    except ZeroDivisionError:
        prec = "NULL"
        
    try:
        fn_num = len([ r for r in grecs if r in intrinsic_crecs ])
    except ZeroDivisionError:
        fn_num = "NULL"
        
    # Stat missing grecs cause
    rough_crecs = [ vcf_rec(r, onlysite=True) for r in raw_crecs if "LIKELY_INTRINSIC" not in r.FILTER ]
    rough_grecs = [ vcf_rec(r, onlysite=True) for r in raw_grecs ]
    tp_crecs = [ r for r in rough_crecs if r in rough_grecs ]
    tp_grecs = [ r for r in rough_grecs if r in rough_crecs ]
    
    tp_map = []
    for grec in tp_grecs:
        same_site_crecs = [ r for r in tp_crecs if r == grec ]
        assert len(same_site_crecs) <= 2
        if len(same_site_crecs) == 1:
            tp_map.append(tuple([grec, same_site_crecs[0]]))
        elif len(same_site_crecs) == 0:
            continue
        else:
            crec1 = same_site_crecs[0]
            crec2 = same_site_crecs[1]
            c1_genotype = crec1.genotype if crec1.genotype != "miss" else np.nan
            c2_genotype = crec2.genotype if crec2.genotype != "miss" else np.nan
            if c1_genotype in [ c2_genotype ]:
                crec = crec1 if "GATK" in crec1.filters else crec2
                tp_map.append(tuple([grec, crec]))
            elif c1_genotype in [1,2] and c2_genotype in [1,2]:
                crec = crec1 if "GATK" in crec1.filters else crec2
                tp_map.append(tuple([grec, crec]))
            else:
                gt_map = { c1_genotype: crec1, c2_genotype: crec2 }
                tp_map.append(tuple([grec, gt_map[max([c1_genotype, c2_genotype])]]))
            
    
    # Stats the vars having matching site but mistakenly called homoalt to a hetero var
    het_to_hom = len([ t for t in tp_map if t[0].genotype == 1 and t[1].genotype == 2])
    gatk_het_to_hom = len([ t for t in tp_map if t[0].genotype == 1 and t[1].genotype == 2 and "GATK" in t[1].filters ])
    # Stats the vars having matching site but mistakenly called hetero to a homoalt var 
    hom_to_het = len([ t for t in tp_map if t[0].genotype == 2 and t[1].genotype == 1])
    gatk_hom_to_het = len([ t for t in tp_map if t[0].genotype == 2 and t[1].genotype == 1 and "GATK" in t[1].filters ])
    # Stats for the tp with the same genotype
    accurate_GT_call = len([ t for t in tp_map if t[1].genotype in [1,2] and t[0].genotype in [t[1].genotype] ])
    gatk_accurate_GT = len([ t for t in tp_map if t[1].genotype in [1,2] and t[0].genotype in [t[1].genotype] and "GATK" in t[1].filters ])
    SDrecall_accurate_GT = len([ t for t in tp_map if t[1].genotype in [1,2] and t[0].genotype in [t[1].genotype] and "SDrecall" in t[1].filters ])
    
    # Stats the vars that not covered
    uncovered = len([ t for t in tp_map if t[0].genotype in [1,2] and t[1].genotype in ["miss"] ] + [r for r in rough_grecs if r not in rough_crecs])
    # Stats the vars covered but ref
    missed = len([ t for t in tp_map if t[0].genotype in [1,2] and t[1].genotype == 0 and t[1].aa_count in range(1,300) ])
    missed_lowaa = len([ t for t in tp_map if t[0].genotype in [1,2] and t[1].genotype == 0 and t[1].aa_count in [0,np.nan] ])
        
    fp_num = len(set([vcf_rec(r.__orirecord__, onlysite=True) for r in crecs if r not in grecs and "SDrecall" in r.filters]))
    
    # Prepare for the note field
    note = ""
    if exonic:
        note = note + "exonic" + "_"
    if gatk:
        note = note + "gatk" + "_"
    if rough:
        note = note + "onlysite" + "_"
    note = note.strip("_")
    
    # Prepare for uncommon varset
    uncommon_crecs = [ r for r in crecs if "COM_or_SE" not in r.filters and "gnomAD_common" not in r.filters ]
    uncommon_crecs_num = len(uncommon_crecs)
    uncommon_grecs = [ r for r in grecs if "COM_or_SE" not in r.filters and "gnomAD_common" not in r.filters ]
    uncommon_grecs_num = len(uncommon_grecs)
    uncommon_fp_num = len([r for r in uncommon_crecs if r not in uncommon_grecs and "SDrecall" in r.filters ])
    try:
        uncommon_recall = len([ r for r in uncommon_grecs if r in uncommon_crecs ])/len(uncommon_grecs)
    except ZeroDivisionError:
        uncommon_recall = "NULL"
    try:
        uncommon_prec = len([r for r in uncommon_crecs if r in uncommon_grecs])/len(uncommon_crecs)
    except ZeroDivisionError:
        uncommon_prec = "NULL"
    
    # Prepare tab creation sql query
    tab_name = "benchmarks"
    tab_create_sql = create_sql_tab(tab_name=tab_name,
                                    keys = ["samp_id", "note"],
                                    samp_id = "TEXT",
                                    note = "TEXT",
                                    called_vcf = "TEXT",
                                    golden_vcf = "TEXT",
                                    crecs_num = "INTEGER",
                                    intrinsic_crecs_num = "INTEGER",
                                    grecs_num = "INTEGER",
                                    fn_num = "INTEGER",
                                    recall = "REAL",
                                    prec = "REAL",
                                    fp_num = "INTEGER",
                                    accurate_GT_call = "INTEGER",
                                    gatk_het2hom = "INTEGER",
                                    gatk_hom2het = "INTEGER", 
                                    gatk_accurate_GT = "INTEGER", 
                                    SDrecall_accurate_GT = "INTEGER", 
                                    het2hom = "INTEGER",
                                    hom2het = "INTEGER",
                                    uncovered = "INTEGER",
                                    missed_with_alt = "INTEGER",
                                    missed_no_alt = "INTEGER",
                                    uncommon_crecs_num = "INTEGER",
                                    uncommon_grecs_num = "INTEGER",
                                    uncommon_recall = "REAL",
                                    uncommon_prec = "REAL",
                                    uncommon_fp_num = "INTEGER")
    
    mod_sql_query = compose_mod_sql_query(tabname=tab_name,
                                          keys = ["samp_id", "note"],
                                          samp_id = samp_id,
                                          note = note,
                                          called_vcf = called_vcf,
                                          golden_vcf = golden_vcf,
                                          crecs_num = crecs_num,
                                          intrinsic_crecs_num = intrinsic_crecs_num,
                                          grecs_num = grecs_num,
                                          fn_num = fn_num,
                                          recall = recall,
                                          prec = prec,
                                          fp_num = fp_num,
                                          accurate_GT_call = accurate_GT_call,
                                          gatk_het2hom = gatk_het_to_hom,
                                          gatk_hom2het = gatk_hom_to_het,
                                          gatk_accurate_GT = gatk_accurate_GT, 
                                          SDrecall_accurate_GT = SDrecall_accurate_GT, 
                                          het2hom = het_to_hom,
                                          hom2het = hom_to_het,
                                          uncovered = uncovered,
                                          missed_with_alt = missed,
                                          missed_no_alt = missed_lowaa,
                                          uncommon_crecs_num = uncommon_crecs_num,
                                          uncommon_grecs_num = uncommon_grecs_num,
                                          uncommon_recall = uncommon_recall,
                                          uncommon_prec = uncommon_prec,
                                          uncommon_fp_num = uncommon_fp_num)
    
    meta_db_path = meta_file[:-4] + ".sqlite3"
    db_recording_benchmark(dbpath = meta_db_path,
                           table_map = {tab_name: tab_create_sql},
                           upsert_sql = [mod_sql_query],
                           output_tab_path = {tab_name: meta_file})
    




class vcf_rec(object):
    '''
    This class is used as a node object in Graph 
    '''
    def __init__(self, vcfpy_record, caller=".", rough = False, onlysite=False):
        self.__chrom__ = str(vcfpy_record.CHROM)
        self.__start__ = int(vcfpy_record.POS)
        try:
            self.__end__ = int(vcfpy_record.INFO['END'])
        except:
            self.__end__ = int(vcfpy_record.POS)
        self.__orirecord__ = vcfpy_record
        self.infos = vcfpy_record.INFO
        self.__ref__ = str(vcfpy_record.REF).upper()
        try:
            self.depth = float(vcfpy_record.calls[0].data["DP"])
        except IndexError:
            self.depth = 0
        except KeyError:
            self.depth = np.nan
        try:
            self.genotype = vcfpy_record.calls[0].gt_type
            if self.depth == 0:
                self.genotype = "miss"
            if self.depth in [np.nan]:
                self.genotype = "miss"
        except IndexError:
            self.genotype = "miss"
        else:
            if "AD" in vcfpy_record.calls[0].data:
                ad_values = vcfpy_record.calls[0].data["AD"]
                if ad_values[0] > 0 and ad_values[1] == 0:
                    self.genotype = 0
            if self.genotype is None:
                self.genotype = "miss"
        try:
            if "AD" in vcfpy_record.calls[0].data:
                self.aa_count = vcfpy_record.calls[0].data["AD"][1]
            elif self.depth > 0:
                self.aa_count = 0
            else:
                self.aa_count = np.nan
        except IndexError:
            self.aa_count = np.nan
            
        try:
            self.__alt__ = str(vcfpy_record.ALT[0]).upper()
        except IndexError as ie:
            # logger.warning("This vcfpy record seems being a ref record. Let's take a look at it: \n{}".format(vcfpy_record))
            self.__alt__ = self.__ref__
        self.__caller__ = str(caller)
        self.__samples__ = list(vcfpy_record.call_for_sample.keys())
        self.filters = vcfpy_record.FILTER
        self.ids = vcfpy_record.ID
        self.rough = rough
        self.onlysite = onlysite

    def __hash__(self):
        if self.onlysite:
            return hash(self.__chrom__) ^ hash(self.__start__) ^ hash(self.__ref__) ^ hash(self.__alt__) ^ hash(self.__caller__) ^ hash(",".join(self.__samples__))
        else:
            return hash(self.__chrom__) ^ hash(self.__start__) ^ hash(self.__ref__) ^ hash(self.__alt__) ^ hash(self.__caller__) ^ hash(",".join(self.__samples__)) ^ hash(self.genotype)

    def __eq__(self, other):
        '''
        Test if this object is equal to other object
        '''
        if self.__chrom__ != other.__chrom__:
            return False
        elif self.onlysite:
            if self.__start__ == other.__start__ and self.__ref__ == other.__ref__ and self.__alt__ == other.__alt__ and self.__samples__ == other.__samples__:
                return True
            else:
                return False
        elif self.rough:
            if self.__start__ == other.__start__ and self.__ref__ == other.__ref__ and self.__alt__ == other.__alt__ and self.__samples__ == other.__samples__ and self.genotype != 0 and other.genotype != 0 and self.genotype != "miss" and other.genotype != "miss":
                return True
            else:
                return False
        else:
            if self.__start__ == other.__start__ and self.__ref__ == other.__ref__ and self.__alt__ == other.__alt__ and self.__samples__ == other.__samples__ and self.genotype == other.genotype and self.genotype != 0 and other.genotype != 0 and self.genotype != "miss" and other.genotype != "miss":
                return True
            else:
                return False
            

    def __str__(self):
        '''
        Convert this object to string
        '''
        return "chr:{},start:{},end:{},ref:{},alt:{},caller:{},gt:{}".format(self.__chrom__, self.__start__, self.__end__, self.__ref__, self.__alt__, self.__caller__, self.genotype)

    def __repr__(self):
        '''
        Representation for printing
        '''
        return "chr:{},start:{},end:{},ref:{},alt:{},caller:{},gt:{}".format(self.__chrom__, self.__start__, self.__end__, self.__ref__, self.__alt__, self.__caller__, self.genotype)
        
        
        
if __name__ == "__main__":
    all_combs = list(it.product([True, False], [True, False], [True, False]))
    for args in all_combs:
        compare_per_pair(sys.argv[1], sys.argv[2], sys.argv[3], rough = args[0], exonic = args[1], gatk = args[2])