#!/usr/bin/env python

# ================================
# @auther: Rongbin Zheng
# @email: Rongbin.Zheng@childrens.harvard.edu
# @date: May 2022
# ================================

import os,sys
import time, re
import pickle as pk
from datetime import datetime
import numpy as np
import pandas as pd
from operator import itemgetter
import scipy
from scipy import sparse
import scanpy as sc
import collections
import multiprocessing
import configparser
import tracemalloc

from matplotlib.backends.backend_pdf import PdfPages

import mebocost.MetEstimator as ME
import mebocost.crosstalk_calculator as CC
import mebocost.crosstalk_plots as CP
import mebocost.pathway_enrichment as PE
import mebocost.pathway_plot as PP


"""
linking input and out 
"""

def info(string):
    """
    print information
    """
    today = datetime.today().strftime("%B %d, %Y")
    now = datetime.now().strftime("%H:%M:%S")
    current_time = today + ' ' + now
    print("[{}]: {}".format(current_time, string))


def _correct_colname_meta_(scRNA_meta, cellgroup_col=[]):
    """
    sometime the column names have different
    """
#     print(scRNA_meta)
    if scRNA_meta is None or scRNA_meta is pd.DataFrame:
        raise KeyError('Please provide cell_ann data frame!')
    
    if cellgroup_col:
        ## check columns names
        for x in cellgroup_col:
            if x not in scRNA_meta.columns.tolist():
                info('ERROR: given cell group identifier {} is not in meta table columns'.format(x))
                raise ValueError('given cell group identifier {} is not in meta table columns'.format(x))
        ## get cell group name
        scRNA_meta['cell_group'] = scRNA_meta[cellgroup_col].astype('str').apply(lambda row: '_'.join(row), axis = 1).tolist()
    else:
        info('no cell group given, try to search cluster and cell_type')
        col_names = scRNA_meta.columns.tolist()
        if 'cell_type' in col_names:
            pass
        elif 'cell_type' not in col_names and 'Cell_Type' in col_names:
            scRNA_meta.columns = ['cell_type' if x.upper() == 'CELL_TYPE' else x for x in col_names]
        elif 'cell_type' not in col_names and 'celltype' in col_names:
            scRNA_meta.columns = ['cell_type' if x.upper() == 'CELLTYPE' else x for x in col_names]
        elif 'cell_type' not in col_names and 'CellType' in col_names:
            scRNA_meta.columns = ['cell_type' if x.upper() == 'CELL TYPE' else x for x in col_names]
        else:
            info('ERROR: "cell_type" not in scRNA meta column names, will try cluster')
            if 'cluster' not in col_names and 'Cluster' in col_names:
                scRNA_meta.columns = ['cluster' if x.upper() == 'CLUSTER' else x for x in col_names]
            else:
                raise KeyError('cluster cannot find in the annotation, and cell_group does not specified')
            raise KeyError('cell_type cannot find in the annotation, and cell_group does not specified'.format(x))
        
        if 'cell_type' in scRNA_meta.columns.tolist():
            scRNA_meta['cell_group'] = scRNA_meta['cell_type'].tolist()
        elif 'cluster' in scRNA_meta.columns.tolist():
            scRNA_meta['cell_group'] = scRNA_meta['cluster'].tolist()
        else:
            raise KeyError('Please a group_col to group single cell')
    return(scRNA_meta)

def _read_config(conf_path):
    """
    read config file
    """
    #read config
    cf = configparser.ConfigParser()
    cf.read(conf_path)
    config = cf._sections
    # remove the annotation:
    for firstLevel in config.keys():
        for secondLevel in config[firstLevel]:
            if '#' in config[firstLevel][secondLevel]:
                config[firstLevel][secondLevel] = config[firstLevel][secondLevel][:config[firstLevel][secondLevel].index('#')-1].rstrip()
    return(config)

def load_obj(path):
    """
    read mebocost object
    """
    f = open(path, 'rb')
    obj_vars = pk.load(f)
    f.close()
    keys = list(obj_vars.keys())
    mebocost_obj = create_obj(exp_mat = obj_vars['exp_mat'],
                        adata = obj_vars['adata'],
                        cell_ann = obj_vars['cell_ann'],
                        group_col = obj_vars['group_col'] if 'group_col' in keys else ['celltype'],
                        config_path = obj_vars['config_path'],
                        met_enzyme = obj_vars['met_enzyme'],
                        met_sensor = obj_vars['met_sensor'],
                        met_ann = obj_vars['met_ann'], 
                        scFEA_ann = obj_vars['scFEA_ann'],
                        compass_met_ann = obj_vars['compass_met_ann'],
                        compass_rxn_ann = obj_vars['compass_rxn_ann'],
                        gene_network = obj_vars['gene_network'],
                        gmt_path = obj_vars['gmt_path']
                       )
    
    mebocost_obj.exp_mat_indexer = obj_vars['exp_mat_indexer'] if 'exp_mat_indexer' in keys else []
    mebocost_obj.exp_mat_columns = obj_vars['exp_mat_columns'] if 'exp_mat_columns' in keys else []
    mebocost_obj.avg_exp = obj_vars['avg_exp'] if 'avg_exp' in keys else pd.DataFrame()
    mebocost_obj.avg_exp_indexer = obj_vars['avg_exp_indexer'] if 'avg_exp_indexer' in keys else []
    mebocost_obj.avg_exp_columns = obj_vars['avg_exp_columns'] if 'avg_exp_columns' in keys else []
    mebocost_obj.avg_exp_columns = obj_vars['avg_exp_columns'] if 'avg_exp_columns' in keys else []
    mebocost_obj.met_mat = obj_vars['met_mat'] if 'met_mat' in keys else pd.DataFrame()
    mebocost_obj.met_mat_indexer = obj_vars['met_mat_indexer'] if 'met_mat_indexer' in keys else pd.DataFrame()
    mebocost_obj.met_mat_columns = obj_vars['met_mat_columns'] if 'met_mat_columns' in keys else pd.DataFrame()
    mebocost_obj.avg_met = obj_vars['avg_met'] if 'avg_met' in keys else pd.DataFrame()
    mebocost_obj.avg_met_indexer = obj_vars['avg_met_indexer'] if 'avg_met_indexer' in keys else pd.DataFrame()
    mebocost_obj.avg_met_columns = obj_vars['avg_met_columns'] if 'avg_met_columns' in keys else pd.DataFrame()

    mebocost_obj.species = obj_vars['species'] if 'species' in keys else 'human'
    mebocost_obj.met_est = obj_vars['met_est'] if 'met_est' in keys else 'mebocost'
    mebocost_obj.met_pred = obj_vars['met_pred'] if 'met_pred' in keys else pd.DataFrame()
    mebocost_obj.cutoff_exp = obj_vars['cutoff_exp'] if 'cutoff_exp' in keys else 0
    mebocost_obj.cutoff_met = obj_vars['cutoff_met'] if 'cutoff_met' in keys else 0
    mebocost_obj.cutoff_prop = obj_vars['cutoff_prop'] if 'cutoff_prop' in keys else 0.1
    mebocost_obj.sensor_type = obj_vars['sensor_type'] if 'sensor_type' in keys else ['Receptor', 'Transporter', 'Nuclear Receptor']
    mebocost_obj.thread = obj_vars['thread'] if 'thread' in keys else 1
    mebocost_obj.commu_res = obj_vars['commu_res'] if 'commu_res' in keys else pd.DataFrame()
    mebocost_obj.original_result = obj_vars['original_result'] if 'original_result' in keys else pd.DataFrame()
    mebocost_obj.commu_bg = obj_vars['commu_bg'] if 'commu_bg' in keys else dict()
    mebocost_obj.exp_prop = obj_vars['exp_prop'] if 'exp_prop' in keys else pd.DataFrame()
    mebocost_obj.met_prop = obj_vars['met_prop'] if 'met_prop' in keys else pd.DataFrame()
    ## check pathway enrich
    mebocost_obj.enrich_result = obj_vars['enrich_result'] if 'enrich_result' in keys else dict() 
    return mebocost_obj


def save_obj(obj, path = 'mebocost_result.pk', filetype = 'pickle'):
    """
    save object to pickle
    """
    if filetype == 'pickle':
        obj_vars = vars(obj)
        out = open(path, 'wb')
        pk.dump(obj_vars, out)
        out.close()
    

class create_obj:
    """
    MEBOCOST for predicting metabolite-based cell-cell communication. The modules of the package include communication inference, communication visualization, pathway inference, pathway visualization.

    Params
    -------
    exp_mat
        python pandas data frame, single cell expression matrix, rows are genes, columns are cells
        'exp_mat' is a exclusive parameter to 'adata'
    adata
        scanpy adata object, the expression will be extracted, 'adata' is an exclusive parameter to 'exp_mat'
    cell_ann
        data frame, cell annotation information, cells in row names
    group_col
        a list, specify the column names in 'cell_ann' for grouping cells, by default 'cell_type' or 'cluster' will be detected and used
    species
        human or mouse, this determines which database will be used in our collection

    met_est
        the method for estimating metabolite level in cell, should be one of:
        mebocost: estimated by the enzyme network related to the metabolite
        scFEA-flux: flux result of published software scFEA (https://pubmed.ncbi.nlm.nih.gov/34301623/)
        scFEA-balance: balance result of published software scFEA (https://pubmed.ncbi.nlm.nih.gov/34301623/)
        compass-reaction: reaction result of published software Compass (https://pubmed.ncbi.nlm.nih.gov/34216539/)
        compass-uptake: uptake result of published software Compass (https://pubmed.ncbi.nlm.nih.gov/34216539/)
        compass-secretion: secretion result of published software Compass (https://pubmed.ncbi.nlm.nih.gov/34216539/)
    met_pred
        data frame, if scFEA or Compass is used to impute the metabolite level in cells, please provide the original result from scFEA or Compass, cells in row names, metabolite/reaction/module in column names, 
        Noted that this parameter will be ignored if 'met_est' was set as mebocost.

    config_path
        str, the path for a config file containing the path of files for metabolite annotation, enzyme, sensor, scFEA annotation, compass annotation. These can also be specified separately by paramters as following:

        if config_path not given, please set:
    met_enzyme
        data frame, metabolite and gene (enzyme) relationships, required columns include HMDB_ID, gene, direction, for instance:
        
        HMDB_ID     gene                                                direction
        HMDB0003375 Cyp2c54[Unknown]; Cyp2c38[Unknown]; Cyp2c50[Un...   substrate
        HMDB0003375 Cyp2c54[Unknown]; Cyp2c38[Unknown]; Cyp2c50[Un...   substrate
        HMDB0003375 Cyp2c54[Unknown]; Cyp2c38[Unknown]; Cyp2c50[Un...   substrate
        HMDB0003450 Cyp2c54[Unknown]; Cyp2c38[Unknown]; Cyp2c50[Un...   product
        HMDB0003948 Tuba8[Unknown]; Ehhadh[Unknown]; Echs1[Enzyme]...   product

    met_sensor
        data frame, metabolite sensor information, each row is a pair of metabolite and sensor, must include columns  HMDB_ID, Gene_name, Annotation, for instance:
        
        HMDB_ID Gene_name   Annotation
        HMDB0006247 Abca1   Transporter
        HMDB0000517 Slc7a1  Transporter
        HMDB0000030 Slc5a6  Transporter
        HMDB0000067 Cd36    Transporter
        
    met_ann:
        data frame, the annotation of metabolite collected from HMDB website, these are basic annotation info including HMDB_ID, Kegg_ID, metabolite, etc

    scFEA_ann
        data frame, module annotation of metabolite flux in scFEA, usually is the file at https://github.com/changwn/scFEA/blob/master/data/Human_M168_information.symbols.csv

    compass_met_ann
        data frame, the metabolite annotation used in Compass software, usually is the file at https://github.com/YosefLab/Compass/blob/master/compass/Resources/Recon2_export/met_md.csv

    compass_rxn_ann
        data frame, the reaction annotation used in Compass software, usually is the file at https://github.com/YosefLab/Compass/blob/master/compass/Resources/Recon2_export/rxn_md.csv

    gene_network
        data frame, gene by gene matrix, the value represent the association between two genes, will be used to evaluate downstream effect of the communication

    gmt_path
        a path, this parameter can be provided in config file and given by config_path. Only set this when you do not pass config_path parameter in. The gmt file contains pathway gene list, will be used in pathway inference module, the details of GMT format could be found at https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#:~:text=The+GMT+file+format+is,genes+in+the+gene+set. 

    cutoff_exp
        float, used to filter out cells which are lowly expressed for the given gene

    cutoff_met
        float, used to filter out cells which are lowly abundant of the given metabolite

    cutoff_prop
        float from 0 to 1, used to filter out metabolite or genes if the proportion of their abundant cells less than the cutoff

    sensor_type
        a list, provide a list of sensor type that will be used in the communication modeling, must be one or more from ['Receptor', 'Transporter', 'Nuclear Receptor'], default is all the three

    thread
        int, number of cores used for running job, default 1
        
    """
    def __init__(self,  
                exp_mat=None, 
                adata=None, 
                cell_ann=None,
                group_col=[],
                species = 'human',

                met_est=None,
                met_pred=pd.DataFrame(), 

                config_path=None,
                met_enzyme=pd.DataFrame(),
                met_sensor=pd.DataFrame(),
                met_ann=pd.DataFrame(), 
                scFEA_ann=pd.DataFrame(),
                compass_met_ann=pd.DataFrame(),
                compass_rxn_ann=pd.DataFrame(),
                gene_network=pd.DataFrame(),
                gmt_path = None,

                cutoff_exp=0,
                cutoff_met=0,
                cutoff_prop=0.25,

                sensor_type=['Receptor', 'Transporter', 'Nuclear Receptor'],
                thread = 1
                ):
        tic = time.time()

        self.exp_mat = exp_mat
        self.adata = adata
        ## check cell group information
        ## add a column "cell_group" if successfull
        if (self.exp_mat is None and cell_ann is None) and (self.adata is not None):
            cell_ann = adata.obs.copy()
        self.group_col = ['cell_type', 'cluster'] if not group_col else group_col
        self.cell_ann = _correct_colname_meta_(cell_ann, cellgroup_col = self.group_col)
        self.species = species

        self.met_est = 'mebocost' if not met_est else met_est # one of [scFEA-flux, scFEA-balance, compass-reaction, compass-uptake, compass-secretion]
        self.met_pred = met_pred

        ## the path of config file
        self.config_path = config_path
        ## genes (enzyme) related to met
        self.met_enzyme = met_enzyme
        ## gene name in metaboltie sensor
        self.met_sensor = met_sensor
        ## met basic ann
        self.met_ann = met_ann
        ## software ann
        self.scFEA_ann = scFEA_ann
        self.compass_met_ann = compass_met_ann
        self.compass_rxn_ann = compass_rxn_ann
        ## gene network
        self.gene_network = gene_network
        if not self.config_path and (self.met_sensor is None or self.met_sensor.shape[0] == 0):
            raise KeyError('Please either provide config_path or a data frame of met_enzyme, met_sensor, met_ann, gene_network, etc')

        ## cutoff for expression, metabolite, and proportion of cells
        self.cutoff_exp = cutoff_exp
        self.cutoff_met = cutoff_met
        self.cutoff_prop = cutoff_prop
        self.sensor_type = sensor_type
        self.gmt_path = gmt_path
        self.thread = thread
        
        ## ============== initial ===========

        if self.exp_mat is None and self.adata is None:
            raise ValueError('ERROR: please provide expression matrix either from exp_mat or adata (scanpy object)')  
        elif self.exp_mat is None and self.adata is not None:
            ## check the adata object
            ngene = len(self.adata.var_names)
            ncell = len(self.adata.obs_names)
            info('We get expression data with {n1} genes and {n2} cells.'.format(n1 = ngene, n2 = ncell))
            if ngene < 5000:
                info('scanpy object contains less than 5000 genes, please make sure you are using raw.to_adata()')
            self.exp_mat = sparse.csc_matrix(self.adata.X.T)
            self.exp_mat_indexer = self.adata.var_names
            self.exp_mat_columns = self.adata.obs_names
            self.adata = None
        else:
            if 'scipy.sparse' in str(type(self.exp_mat)):
                ## since the scipy version problem leads to the failure of using sparse.issparse
                ## use a simple way to check!!!
                #sparse.issparse(self.exp_mat):
                pass 
            elif type(self.exp_mat) is type(pd.DataFrame()):
                self.exp_mat_indexer = self.exp_mat.index ## genes
                self.exp_mat_columns = self.exp_mat.columns ## columns
                self.exp_mat = sparse.csc_matrix(self.exp_mat)
                ngene, ncell = self.exp_mat.shape
                info('We get expression data with {n1} genes and {n2} cells.'.format(n1 = ngene, n2 = ncell))
            else:
                info('ERROR: cannot read the expression matrix, please provide pandas dataframe or scanpy adata')

        ## end preparation
        toc = time.time()
        info('Data Preparation Done in {:.4f} seconds'.format(toc-tic))
    
    def _load_config_(self):
        """
        load config and read data from the given path based on given species
        """
        ## the path of config file
        info('Load config and read data based on given species [%s].'%(self.species))
        if self.config_path:
            if not os.path.exists(self.config_path):
                raise KeyError('ERROR: the config path is not exist!')
            config = _read_config(conf_path = self.config_path)
            ## common
            self.met_ann = pd.read_csv(config['common']['hmdb_info_path'], sep = '\t')
            if self.met_est.startswith('scFEA'):
                    self.scFEA_ann = pd.read_csv(config['common']['scfea_info_path'], index_col = 0)
            if self.met_est.startswith('compass'):
                self.compass_met_ann = pd.read_csv(config['common']['compass_met_ann_path'])
                self.compass_rxn_ann = pd.read_csv(config['common']['compass_rxt_ann_path'])
            ## depends on species
            if self.species == 'human':
                self.met_enzyme = pd.read_csv(config['human']['met_enzyme_path'], sep = '\t')
                met_sensor = pd.read_csv(config['human']['met_sensor_path'], sep = '\t')
#                 met_sensor['gene'] = met_sensor['Gene_name'].apply(lambda x: x.split('[')[0])
                self.met_sensor = met_sensor
                self.gene_network = pd.read_csv(config['human']['gene_network_path'], index_col = 0, compression = 'gzip')
                self.gmt_path = config['human']['kegg_gmt_path']
            elif self.species == 'mouse':
                self.met_enzyme = pd.read_csv(config['mouse']['met_enzyme_path'], sep = '\t')
                met_sensor = pd.read_csv(config['mouse']['met_sensor_path'], sep = '\t')
#                 met_sensor['gene'] = met_sensor['Gene_name'].apply(lambda x: x.split('[')[0])
                self.met_sensor = met_sensor
                self.gene_network = pd.read_csv(config['mouse']['gene_network_path'], index_col = 0, compression = 'gzip')
                self.gmt_path = config['mouse']['kegg_gmt_path']
            else:
                raise KeyError('Species should be either human or mouse!')
                            ## check row and columns, we expect rows are genes, columns are cells
                if len(set(self.met_sensor['Gene_name'].tolist()) & set(self.exp_mat_indexer.tolist())) < 10 and len(set(self.met_sensor['Gene_name'].tolist()) & set(self.exp_mat_columns.tolist())) < 10:
                    raise KeyError('it looks like that both the row and columns are not matching to gene name very well, please check the provided matrix or species!')
                if len(set(self.met_sensor['Gene_name'].tolist()) & set(self.exp_mat_indexer.tolist())) < 10 and len(set(self.met_sensor['Gene_name'].tolist()) & set(self.exp_mat_columns.tolist())) > 10:
                    info('it is likely the columns of the exp_mat are genes, will transpose the matrix')
                    self.exp_mat = self.exp_mat.T
                    columns = self.exp_mat_indexer.copy()
                    index = self.exp_mat_columns.copy()
                    self.exp_mat_indexer = index
                    self.exp_mat_columns = columns
        else:
            info('please provide config path')

    def estimator(self):
        """
        estimate of metabolite level in cells using the expression of related enzymes
        """
        info('Estimtate metabolite level using %s'%self.met_est)
        mtd = self.met_est

        if mtd == 'mebocost':
            met_mat, met_indexer, met_columns = ME._met_from_enzyme_est_(exp_mat=self.exp_mat, 
                                                   indexer = self.exp_mat_indexer,
                                                   columns = self.exp_mat_columns,
                                                    met_gene=self.met_enzyme, 
                                                    method = 'mean')
        elif mtd == 'scFEA-flux':
            met_mat = ME._scFEA_flux_est_(scFEA_pred = self.met_pred, 
                                            scFEA_info=self.scFEA_ann, 
                                            hmdb_info=self.met_ann)
        elif mtd == 'scFEA-balance':
            met_mat = ME._scFEA_balance_est_(scFEA_pred = self.met_pred, 
                                                scFEA_info=self.scFEA_ann, 
                                                hmdb_info=self.met_ann)
        elif mtd == 'compass-reaction':
            met_mat = ME._compass_react_est_(compass_pred=self.met_pred, 
                                                compass_react_ann=self.compass_rxn_ann, 
                                                compass_met_ann=self.compass_met_ann, 
                                                hmdb_info=self.met_ann)
        else:
            raise KeyError('Please specify "met_est" to be one of [mebocost, scFEA-flux, scFEA-balance, compass-reaction, compass-uptake, compass-secretion]')
        
        self.met_mat = sparse.csc_matrix(met_mat)
        self.met_mat_indexer = np.array(met_indexer)
        self.met_mat_columns = np.array(met_columns)
#         return met_mat


    def infer(self, met_mat=pd.DataFrame(), n_shuffle = 1000, seed = 12345, thread = None):
        """
        excute communication prediction
        met_mat
            data frame, columns are cells and rows are metabolites
        """
        info('Infer communications')
        if met_mat.shape[0] != 0: ## if given met_mat in addition
            self.met_mat_indexer = np.array(met_mat.index)
            self.met_mat_columns = np.array(met_mat.columns)
            self.met_mat = sparse.csc_matrix(met_mat)
        ## focus on met and gene of those are in the data matrix
        met_sensor = self.met_sensor[self.met_sensor['Gene_name'].isin(self.exp_mat_indexer) & 
                                     self.met_sensor['HMDB_ID'].isin(self.met_mat_indexer)]
        self.met_sensor = met_sensor

        ## init
        cobj = CC.InferComm(exp_mat = self.exp_mat,
                            exp_mat_indexer = self.exp_mat_indexer, 
                            exp_mat_columns = self.exp_mat_columns,
                            avg_exp = self.avg_exp,
                            avg_exp_indexer = self.avg_exp_indexer,
                            avg_exp_columns = self.avg_exp_columns,
                            met_mat = self.met_mat,
                            met_mat_indexer = self.met_mat_indexer,
                            met_mat_columns = self.met_mat_columns,
                            avg_met = self.avg_met,
                            avg_met_indexer = self.avg_met_indexer,
                            avg_met_columns = self.avg_met_columns,
                            cell_ann = self.cell_ann,
                            met_sensor = self.met_sensor,
                            sensor_type = self.sensor_type,
                            thread = thread
                           )

        commu_res_df, commu_res_bg = cobj.pred(n_shuffle = n_shuffle, seed = seed)
    
        ## add metabolite name
        hmdbid_to_met = {}
        for Id, met in self.met_ann[['HMDB_ID', 'metabolite']].values.tolist():
            hmdbid_to_met[Id] = met
        ## add name
        commu_res_df['Metabolite_Name'] = list(map(lambda x: hmdbid_to_met.get(x) if x in hmdbid_to_met else None,
                                                   commu_res_df['Metabolite']))

        ## add annotation
        sensor_to_ann = {}
        for s, a in self.met_sensor[['Gene_name', 'Annotation']].values.tolist():
            sensor_to_ann[s] = a
        commu_res_df['Annotation'] = list(map(lambda x: sensor_to_ann.get(x) if x in sensor_to_ann else None,
                                              commu_res_df['Sensor']))
        
        return commu_res_df, commu_res_bg


    def _filter_lowly_aboundant_(self, 
                                 pvalue_res,
                                 cutoff_prop,
                                 met_prop=None,
                                 exp_prop=None,
                                 min_cell_number=50
                                ):
        """
        change p value to 1 if either metabolite_prop or transporter_prop equal to 0 
        (meaning that no metabolite or transporter level in the cluster)
        """
        ## add the metabolite abudance proportion
        if met_prop is not None:
            pvalue_res['metabolite_prop_in_sender'] = [met_prop.loc[s, m] for s, m in pvalue_res[['Sender', 'Metabolite']].values.tolist()]
        ## add the metabolite abudance proportion
        if exp_prop is not None:
            pvalue_res['sensor_prop_in_receiver'] = [exp_prop.loc[r, s] for r, s in pvalue_res[['Receiver', 'Sensor']].values.tolist()]
        
        if 'original_result' not in list(vars(self)):
            self.original_result = pvalue_res.copy()
        ## minimum cell number
        cell_count = pd.Series(dict(collections.Counter(self.cell_ann['cell_group'].tolist())))
        bad_cellgroup = cell_count[cell_count<min_cell_number].index.tolist() 
        
        info('Set p value and fdr to 1 if sensor or metaboltie expressed cell proportion less than {}'.format(cutoff_prop))
        bad_index = np.where((pvalue_res['metabolite_prop_in_sender'] <= cutoff_prop) |
                             (pvalue_res['sensor_prop_in_receiver'] <= cutoff_prop) |
                             (pvalue_res['Commu_Score'] < 0) |
                             (pvalue_res['Sender'].isin(bad_cellgroup)) | 
                             (pvalue_res['Receiver'].isin(bad_cellgroup))
                            )[0]
        if len(bad_index) > 0:
            pval_index = np.where(pvalue_res.columns.str.endswith('_pval'))[0]
            pvalue_res.iloc[bad_index, pval_index] = 1 # change to 1
            fdr_index = np.where(pvalue_res.columns.str.endswith('_fdr'))[0]
            pvalue_res.iloc[bad_index, fdr_index] = 1 # change to 1
        
        ## norm communication score
        pvalue_res['Commu_Score'] = pvalue_res['Commu_Score']/np.array(pvalue_res['bg_mean']).clip(min = 0.05)
        
        ## reorder columns
        columns = ['Sender', 'Metabolite', 'Metabolite_Name', 
                   'Receiver', 'Sensor', 'Commu_Score', 
                   'metabolite_prop_in_sender',
                   'sensor_prop_in_receiver', 
#                    'bg_mean', 'bg_std',
                   'ztest_stat', 'ztest_pval', 'ttest_stat',
                   'ttest_pval', 'ranksum_test_stat', 'ranksum_test_pval',
                   'permutation_test_stat', 'permutation_test_pval',
                   'ztest_fdr', 'ttest_fdr', 'ranksum_test_fdr',
                   'permutation_test_fdr']
        get_columns = [x for x in columns if x in pvalue_res.columns.tolist()]
        pvalue_res = pvalue_res.reindex(columns = get_columns).sort_values('permutation_test_fdr')
        
        return(pvalue_res)


    def _check_aboundance_(self, cutoff_exp=None, cutoff_met=None):
        """
        check the aboundance of metabolite or transporter expression in cell clusters,
        return the percentage of cells that meet the given cutoff
        by default, cutoff for metabolite aboundance is 0, expression of transporter is 0
        """
        info('Calculating aboundance of metabolite and sensor expression in cell groups')
        ## this will re-write the begin values
        if not cutoff_exp:
            cutoff_exp = self.cutoff_exp if isinstance(self.cutoff_exp, float) else 0
        if not cutoff_met:
            cutoff_met = self.cutoff_met if isinstance(self.cutoff_met, float) else 0
        ## expression for all transporters
        sensors = self.met_sensor['Gene_name'].unique().tolist()
        info('cutoff_exp: {}'.format(cutoff_exp))
        
        sensor_loc = {g:i for i,g in enumerate(self.exp_mat_indexer) if g in sensors}
        exp_prop = {}
        for x in self.cell_ann['cell_group'].unique().tolist():
            cells = self.cell_ann[self.cell_ann['cell_group'] == x].index.tolist()
            cell_loc = [i for i, c in enumerate(self.exp_mat_columns) if c in cells]
            s = self.exp_mat[list(sensor_loc.values()),:][:,cell_loc]
            exp_prop[x] = pd.Series([v[v>cutoff_exp].shape[1] / v.shape[1] for v in s],
                                   index = list(sensor_loc.keys()))
        exp_prop = pd.DataFrame.from_dict(exp_prop, orient = 'index')
            
#         emat = pd.merge(exp_df.reindex(index = sensors).T, self.cell_ann[['cell_group']],left_index = True, right_index = True)
#         exp_prop = emat.groupby('cell_group').apply(lambda df: df.apply(lambda col: len(col[col>cutoff_exp]) / len(col))) ## the proportion of gene expressing cells
        
        # ====================== #
        info('cutoff_metabolite: {}'.format(cutoff_met))
        ## metabolite aboundance
        metabolites = self.met_sensor['HMDB_ID'].unique().tolist()
        met_prop = {}
        for x in self.cell_ann['cell_group'].unique().tolist():
            cells = self.cell_ann[self.cell_ann['cell_group'] == x].index.tolist()
            cell_loc = [i for i, c in enumerate(self.met_mat_columns) if c in cells]
            m = self.met_mat[:,cell_loc]
            met_prop[x] = pd.Series([v[v>cutoff_met].shape[1] / v.shape[1] for v in m],
                                   index = self.met_mat_indexer.tolist())
        met_prop = pd.DataFrame.from_dict(met_prop, orient = 'index')

        return exp_prop, met_prop ## cell_group x sensor gene, cell_group x metabolite
    
    def _get_gene_exp_(self):
        """
        only sensor and enzyme gene expression are needed for each cells
        """
        sensors = self.met_sensor['Gene_name'].unique().tolist()
        enzymes = []
        for x in self.met_enzyme['gene'].tolist():
            enzymes.extend([i.split('[')[0] for i in x.split('; ')])
        genes = list(set(sensors+enzymes))
        ## gene loc
        gene_loc = np.where(pd.Series(self.exp_mat_indexer).isin(genes))[0]
        
        gene_dat = self.exp_mat[gene_loc].copy()
        ## update the exp_mat and indexer
        self.exp_mat = sparse.csr_matrix(gene_dat)
        self.exp_mat_indexer = self.exp_mat_indexer[gene_loc]
                                   
    def _avg_by_group_(self):
        ## avg exp by cell_group for met sensor
        group_names = self.cell_ann['cell_group'].unique().tolist()
        avg_exp = np.empty(shape = (self.exp_mat.shape[0],0)) ## save exp data

        for x in group_names:
            cells = self.cell_ann[self.cell_ann['cell_group'] == x].index.tolist()
            cell_loc = np.where(pd.Series(self.exp_mat_columns).isin(cells))[0]
            # arithmatic mean
            avg_exp = np.concatenate((avg_exp, self.exp_mat[:,cell_loc].mean(axis = 1)), axis = 1)
        
        self.avg_exp = sparse.csr_matrix(avg_exp)
        self.avg_exp_indexer = np.array(self.exp_mat_indexer)
        self.avg_exp_columns = np.array(group_names)
    
    
    def _avg_met_group_(self):
        """
        take average of sensor expression and metabolite by cell groups
        """
        ## avg met by cell_group for met
        avg_met = np.empty(shape = (self.met_mat.shape[0],0)) ## save exp data
        group_names = self.cell_ann['cell_group'].unique().tolist()

        for x in group_names:
            cells = self.cell_ann[self.cell_ann['cell_group'] == x].index.tolist()
            cell_loc = cell_loc = np.where(pd.Series(self.met_mat_columns).isin(cells))[0]
            ## mean
            avg_met = np.concatenate((avg_met, self.met_mat[:,cell_loc].mean(axis = 1)), axis = 1)

        self.avg_met = sparse.csr_matrix(avg_met)
        self.avg_met_indexer = np.array(self.met_mat_indexer)
        self.avg_met_columns = group_names

        
    def infer_commu(self, 
                      n_shuffle = 1000,
                      seed = 12345, 
                      Return = True, 
                      thread = None,
                      save_permuation = False,
                      min_cell_number = 50
                     ):
        """
        execute mebocost to infer communications

        Params
        -----
        n_shuffle
            int, number of cell label shuffling for generating null distribution when calculating p-value
            
        seed
            int, a random seed for shuffling cell labels, set seed to get reproducable shuffling result 
            
        Return
            True or False, set True to return the communication event in a data frame
            
        thread
            int, the number of cores used in the computing, default None, thread set when create the object has the highest priority to be considered, so only set thread here if you want to make a change
            
        save_permuation
            True or False, set True to save the communication score for each permutation, this could occupy a higher amount of space when saving out, so default is False

        min_cell_number
            int, the cell groups will be excluded and p-value will be replaced to 1 if there are not enough number of cells (less than min_cell_number), default is 50

        """
        tic = time.time()
        today = datetime.today().strftime("%B %d, %Y")
        now = datetime.now().strftime("%H:%M:%S")
        current_time = today + ' ' + now
        self.commu_time_stamp = current_time
        
        tracemalloc.start()
        ## load config
        self._load_config_()
        
        ## take average by cell group, this must be done before extract sensor and enzyme gene expression of cells
        self._avg_by_group_()
        
        ## extract exp data for sensor and enzyme genes for all cells
        self._get_gene_exp_()
        
        ## estimate metabolite
        self.estimator()
        
        ## avg met mat
        self._avg_met_group_()
    
        # running communication inference
        commu_res_df, commu_res_bg = self.infer(
                                                n_shuffle = n_shuffle, 
                                                seed = seed,
                                                thread = self.thread if thread is None else thread ## allow to set thread in this function
                                                )
        ## update self
        self.commu_res = commu_res_df
        if save_permuation:
            self.commu_bg = commu_res_bg
        
        ## check cell proportion
        exp_prop, met_prop = self._check_aboundance_()
        ## update self
        self.exp_prop = exp_prop
        self.met_prop = met_prop 
        
        ## check low and set p val to 1
        commu_res_df_updated = self._filter_lowly_aboundant_(pvalue_res = commu_res_df,
                                                             cutoff_prop = self.cutoff_prop,
                                                             met_prop=self.met_prop, 
                                                             exp_prop=self.exp_prop,
                                                             min_cell_number = min_cell_number)
        ## update self
        self.commu_res = commu_res_df_updated
        
        current, peak = tracemalloc.get_traced_memory()
        
        # stopping the library
        tracemalloc.stop()
        
        toc = time.time()
        info('Prediction Done in {:.4f} seconds'.format(toc-tic))
        info('Memory Usage in Peak {:.2f} GB'.format(peak / 1024 / 1024 / 1024))
        if Return:
            return(commu_res_df_updated)
        
## ============================== pathway inference ============================
    def infer_pathway(self, 
                     pval_method='permutation_test_fdr', 
                     pval_cutoff = 0.05, 
                     commu_score_cutoff = 0, 
                     commu_score_column = 'Commu_Score', 
                     min_term = 15, 
                     max_term = 500,
                     thread = None,
                     sender_focus = [],
                     metabolite_focus = [],
                     sensor_focus = [],
                     receiver_focus = [],
                     Return_res = False, 
                    ):
        """
        execute MEBOCOST to infer communication associated pathways
        
        Param
        ------
        pval_method
            should be one of ['zztest_pval', 'ttest_pval', 'ranksum_test_pval', 'permutation_test_pval', 'zztest_fdr', 'ttest_fdr', 'ranksum_test_fdr', 'permutation_test_fdr'], default is permutation_test_fdr
        
        pval_cutoff
            float, a value in range between 0 and 1, pvalue less than the cutoff considered as significant event
            
        commu_score_cutoff
            float, communication score greater than the cutoff considered as a good event
            
        commu_score_column
            str, a column name in commu_res table, the column will be considered as communication score column and the communication score greater than the cutoff set by commu_score_cutoff considered as a good event

        min_term
            int, the pathway will be included when the number of genes in the pathway greater than min_term, default is 15
            
        max_term
            int, the pathway will be included when the number of genes in the pathway less than max_term, default is 500

        
        thread 
            int, the number of cores used in the computing, default None, thread set when create the object has the highest priority to be considered, so only set thread here if you want to make a change
        
        sender_focus
            a list of sender cell type or cell groups that will be focused in the analysis
            
        metabolite_focus
            a list of metabolite name that will be focused in the analysis
            
        sensor_focus
            a list of sensor name that will be focused in the analysis
            
        receiver_focus
            a list of receiver cell type or cell groups that will be focused in the analysis
        
        Return_res
            True or False, set True to return the pathway enrichment result after running this function. set False to indicate that do not return in this function but the result will be saved in MEBOCOST object, default is False
            
        """
        tic = time.time()
        today = datetime.today().strftime("%B %d, %Y")
        now = datetime.now().strftime("%H:%M:%S")
        current_time = today + ' ' + now
        self.pathway_time_stamp = current_time
        
        ## check communication 
        if 'commu_res' not in dir(self):
            info('Communication result cannot be found, please excute _predict_() function first!')
            return
        
        ## good communications 
        good_commu = self.commu_res[(self.commu_res[pval_method] < pval_cutoff) &
                                    (self.commu_res[commu_score_column] > commu_score_cutoff)]
        ## start a object                      
        eobj = PE.PathwayEnrich(commu_res = good_commu,
                                gene_network=self.gene_network,
                                avg_exp = self.avg_exp,
                                avg_exp_indexer = self.avg_exp_indexer,
                                avg_exp_columns = self.avg_exp_columns,
                                cell_ann = self.cell_ann,
                                gmt_path = self.gmt_path,
                                min_term = min_term,
                                max_term = max_term,
                                thread = self.thread if thread is None else thread
                                )
        
        self.enrich_result = eobj._pred_(pval_method = pval_method, 
                            sensor_in_receiver = True, 
                            sender_to_receiver = True,
                            sender_focus = sender_focus,
                            metabolite_focus = metabolite_focus,
                            sensor_focus = sensor_focus,
                            receiver_focus = receiver_focus,
                            Return = True
                           )
        
#         ## add
#         self.enrich_result = vars(eobj)
        
        toc  = time.time()
        info('Pathway Inference Done in {:.4f} seconds'.format(toc-tic))
        if Return_res:
            return self.enrich_result['sensor_res'], self.enrich_result['cellpair_res']
        

## ============================== communication plot functions ============================
    def eventnum_bar(self,
                    sender_focus = [],
                    metabolite_focus = [],
                    sensor_focus = [],
                    receiver_focus = [],
                    and_or = 'and',
                    pval_method = 'permutation_test_fdr',
                    pval_cutoff = 0.05,
                    comm_score_col = 'Commu_Score',
                    comm_score_cutoff = None,
                    cutoff_prop = None,
                    figsize = 'auto',
                    save = None,
                    show_plot = True,
                    include = ['sender-receiver', 'sensor', 'metabolite', 'metabolite-sensor'],
                    group_by_cell = True,
                    colorcmap = 'tab20',
                    return_fig = False
                  ):
        """
        this function summarize the number of communication events
        
        Params
        ------
        sender_focus
            a list, set a list of sender cells to be focused, only plot related communications
        metabolite_focus
            a list, set a list of metabolites to be focused, only plot related communications
        sensor_focus
            a list, set a list of sensors to be focused, only plot related communications
        receiver_focus
            a list, set a list of receiver cells to be focused, only plot related communications
        and_or
            eithor 'and' or 'or', 'and' for finding communications that meet to all focus, 'or' for union
        pval_method
            should be one of ['zztest_pval', 'ttest_pval', 'ranksum_test_pval', 'permutation_test_pval', 'zztest_fdr', 'ttest_fdr', 'ranksum_test_fdr', 'permutation_test_fdr'], default is permutation_test_fdr
        pval_cutoff
            float, set to filter out non-significant communication events
        figsize
            auto or a tuple of float such as (5.5, 4.2), defualt will be automatically estimate
        save
            str, the file name to save the figure
        show_plot
             True or False, whether print the figure on the screen
        comm_score_col
            column name of communication score, can be Commu_Score
        comm_score_cutoff
            a float, set a cutoff so only communications with score greater than the cutoff will be focused
        cutoff_prop
            a float between 0 and 1, set a cutoff to further filter out lowly abundant cell populations by the fraction of cells expressed sensor genes or metabolite, Note that this parameter will lost the function if cutoff_prop was set lower than the one user set at begaining of running mebocost.infer_commu or preparing mebocost object. This parameter were designed to further strengthen the filtering.
        include
            a list, contains one or more elements from ['sender-receiver', 'sensor', 'metabolite', 'metabolite-sensor'], we try to summarize the number of communications grouping by the given elements, if return_fig set to be True, only provide one for each run.
        group_by_cell
            True or False, only effective for metabolite and sensor summary, True to further label number of communications in cell groups, False to do not do that
        colormap
            only effective when group_by_cell is True, should be a python camp str, default will be 'tab20', or can be a dict where keys are cell group, values are RGB readable color
        return_fig:
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself.
            
        """
        
#         if show_plot is None and self.show_plot is not None:
#             show_plot = self.show_plot
        
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None
            
        fig = CP._eventnum_bar_(commu_res = self.commu_res,
                    sender_focus = sender_focus,
                    metabolite_focus = metabolite_focus,
                    sensor_focus = sensor_focus,
                    receiver_focus = receiver_focus,
                    and_or = and_or,
                    pval_method = pval_method,
                    pval_cutoff = pval_cutoff,
                    comm_score_col = comm_score_col,
                    comm_score_cutoff = comm_score_cutoff,
                    cutoff_prop = cutoff_prop,
                    figsize = figsize,
                    pdf = Pdf,
                    show_plot = show_plot,
                    include = include,
                    group_by_cell = group_by_cell,
                    colorcmap = colorcmap,
                    return_fig = return_fig
                  )
        if save is not None and save is not False:
            Pdf.close()
        if return_fig:
            return(fig)
    
    def histogram(self, 
                    met_name=None,
                    sensor=None,
                    title = '', 
                    save = None, 
                    show_plot = True,
                    bins = 100,
                    alpha = .6,
                    bg_color = 'grey',
                    obs_color = 'red',
                    figsize = (5.5, 4.2),
                    comm_score_col = 'Commu_Score',
                   return_fig = False):
        """
        histogram plot to show the communication score distribution in background and given pairs of metabolite and sensor

        Params
        -----
        met
            str, metabolite name in the communication table
        sensor
            str, gene name of metabolite sensor in communication table
        title
            str, the title for the figure
        save
            str, the file name to save the figure
        show_plot
            True or False, whether print the figure on the screen
        bins
            int, how many bins plot in the histogram, 100 for default
        alpha
            float, set for color transparent, defaut is 0.6
        bg_color
            color for bars of backgroud communication scores
        obs_color
            color for bars of observed communication score
        comm_score_col
            column name of communication score, can be Commu_Score
        return_fig:
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself

        """
#         if show_plot is None and self.show_plot is not None:
#             show_plot = self.show_plot

        commu_mat = self.original_result[(self.original_result['Metabolite_Name'] == met_name) &
                                    (self.original_result['Sensor'] == sensor)]
        if commu_mat.shape[0] == 0:
            info('ERROR: no data found, please check met_name and sensor') 
            return
        ## find met HMDB ID
        hmdbId = commu_mat['Metabolite'].tolist()[0]
        commu_bg = self.commu_bg[hmdbId+'~'+sensor]
        ## create pdf if true
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None

        fig = CP._histogram_(commu_mat=commu_mat, commu_bg=commu_bg, title = title, pdf = Pdf, 
                 show_plot = show_plot, bins = bins, alpha = alpha, bg_color = bg_color,
                 obs_color = obs_color, figsize = figsize, comm_score_col = comm_score_col,
                            return_fig = return_fig)

        if save is not None and save is not False:
            Pdf.close()
        if return_fig:
            return(fig)

    def commu_dotmap(self,
                sender_focus = [],
                metabolite_focus = [],
                sensor_focus = [],
                receiver_focus = [],
                and_or = 'and',
                pval_method='permutation_test_fdr',
                pval_cutoff=0.05, 
                figsize = 'auto',
                cmap = 'Reds',
                node_size_norm = (10, 150),
                save = None, 
                show_plot = True,
                comm_score_col = 'Commu_Score',
                comm_score_cutoff = None,
                cutoff_prop = None,
                swap_axis = False,
                return_fig = False):
        """
        commu_dotmap to show all significant communication events
        
        Params
        -----
        sender_focus
            a list, set a list of sender cells to be focused, only plot related communications
        metabolite_focus
            a list, set a list of metabolites to be focused, only plot related communications
        sensor_focus
            a list, set a list of sensors to be focused, only plot related communications
        receiver_focus
            a list, set a list of receiver cells to be focused, only plot related communications
        and_or
            eithor 'and' or 'or', 'and' for finding communications that meet to all focus, 'or' for union
        pval_method
            should be one of ['zztest_pval', 'ttest_pval', 'ranksum_test_pval', 'permutation_test_pval', 'zztest_fdr', 'ttest_fdr', 'ranksum_test_fdr', 'permutation_test_fdr'], default is permutation_test_fdr
        pval_cutoff
            float, set to filter out non-significant communication events
        figsize
            auto or a tuple of float such as (5.5, 4.2), defualt will be automatically estimate
        cmap
            colormap for dot color, default is Reds
        node_size_norm
            two values in a tuple, used to normalize the dot size, such as (10, 150)
        save
            str, the file name to save the figure
        show_plot
             True or False, whether print the figure on the screen
        comm_score_col
            column name of communication score, can be Commu_Score
        comm_score_cutoff
            a float, set a cutoff so only communications with score greater than the cutoff will be focused
        cutoff_prop
            a float between 0 and 1, set a cutoff to filter out lowly abundant cell populations by the fraction of cells expressed sensor genes or metabolite, Note that this parameter will lost the function if cutoff_prop was set lower than the one user set at begaining of running mebocost.infer_commu or preparing mebocost object. This parameter were designed to further strengthen the filtering.
        return_fig:
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself
        """
        
#         if show_plot is None and self.show_plot is not None:
#             show_plot = self.show_plot

        comm_res = self.commu_res
        ## pdf
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None
        
        fig = CP._commu_dotmap_(comm_res=comm_res, 
                     sender_focus = sender_focus,
                     metabolite_focus = metabolite_focus,
                     sensor_focus = sensor_focus,
                     receiver_focus = receiver_focus,
                     and_or = and_or,
                     pval_method=pval_method, 
                     pval_cutoff=pval_cutoff, 
                     figsize = figsize, 
                     comm_score_col = comm_score_col,
                     comm_score_cutoff = comm_score_cutoff,
                     cutoff_prop = cutoff_prop,
                     cmap = cmap,
                     node_size_norm = node_size_norm,
                     pdf = Pdf, 
                     show_plot = show_plot,
                     swap_axis = swap_axis,
                     return_fig = return_fig
                    )
        if save is not None and save is not False:
            Pdf.close()
        if return_fig:
            return(fig)
        
    def FlowPlot(self, 
                pval_method='permutation_test_fdr',
                pval_cutoff=0.05,
                sender_focus = [],
                metabolite_focus = [],
                sensor_focus = [],
                receiver_focus = [],
                remove_unrelevant = False,
                and_or = 'and',
                node_label_size = 8,
                node_alpha = .8,
                figsize = 'auto',
                node_cmap = 'Set1',
                line_cmap = 'spring_r',
                line_vmin = None,
                line_vmax = None,
                linewidth_norm = (0.1, 1),
                node_size_norm = (10, 150),
                save=None, 
                show_plot = False,
                comm_score_col = 'Commu_Score',
                comm_score_cutoff = None,
                cutoff_prop = None,
                text_outline = False,
                return_fig = False):
        """
        Flow plot to show the communication connections from sender to metabolite, to sensor, to receiver

        Params
        ------
        pval_method
            should be one of ['zztest_pval', 'ttest_pval', 'ranksum_test_pval', 'permutation_test_pval', 'zztest_fdr', 'ttest_fdr', 'ranksum_test_fdr', 'permutation_test_fdr'], default is permutation_test_fdr
        pval_cutoff
            float, set to filter out non-significant communication events
        sender_focus
            a list, set a list of sender cells to be focused, only plot related communications
        metabolite_focus
            a list, set a list of metabolites to be focused, only plot related communications
        sensor_focus
            a list, set a list of sensors to be focused, only plot related communications
        receiver_focus
            a list, set a list of receiver cells to be focused, only plot related communications
        remove_unrelevant
            True or False, set True to hide unrelated nodes 
        and_or
            eithor 'and' or 'or', 'and' for finding communications that meet to all focus, 'or' for union
        node_label_size
            float, font size of text label on node, default will be 8
        node_alpha
            float, set to transparent node color
        figsize
            auto or a tuple of float such as (5.5, 4.2), defualt will be automatically estimate
        node_cmap
            node color map or a four-element list, used to color sender, metabolite, sensor, receiver, set one from https://matplotlib.org/stable/tutorials/colors/colormaps.html
        line_cmap
            line color map, used to indicate the communication score, set one from https://matplotlib.org/stable/tutorials/colors/colormaps.html
        node_size_norm
            two values in a tuple, used to normalize the dot size, such as (10, 150)
        linewidth_norm
            two values in a tuple, used to normalize the line width, such as (0.1, 1)
        save
            str, the file name to save the figure
        show_plot
            True or False, whether print the figure on the screen
        comm_score_col
            column name of communication score, can be Commu_Score
        comm_score_cutoff
            a float, set a cutoff so only communications with score greater than the cutoff will be focused
        cutoff_prop
            a float between 0 and 1, set a cutoff to filter out lowly abundant cell populations by the fraction of cells expressed sensor genes or metabolite, Note that this parameter will lost the function if cutoff_prop was set lower than the one user set at begaining of running mebocost.infer_commu or preparing mebocost object. This parameter were designed to further strengthen the filtering.
        return_fig:
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself
        """
        
#         if show_plot is None and self.show_plot is not None:
#             show_plot = self.show_plot

        comm_res = self.commu_res
        ## pdf
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None

        fig = CP._FlowPlot_(comm_res=comm_res, pval_method=pval_method, pval_cutoff=pval_cutoff, 
                      sender_focus = sender_focus, metabolite_focus = metabolite_focus,
                      sensor_focus = sensor_focus, receiver_focus = receiver_focus, 
                      remove_unrelevant = remove_unrelevant, and_or = and_or,
                      node_label_size = node_label_size, node_alpha = node_alpha, figsize = figsize, 
                      node_cmap = node_cmap, line_cmap = line_cmap, line_vmin = line_vmin,
                      line_vmax = line_vmax, linewidth_norm = linewidth_norm, 
                      node_size_norm = node_size_norm, pdf=Pdf, show_plot = show_plot, 
                      comm_score_col = comm_score_col, comm_score_cutoff = comm_score_cutoff, cutoff_prop = cutoff_prop,
                      text_outline = text_outline, return_fig = return_fig)
        if save is not None and save is not False:
            Pdf.close()
        if return_fig:
            return(fig)
    def count_dot_plot(self, 
                    pval_method='permutation_test_pval', 
                    pval_cutoff=0.05, 
                    cmap='RdBu_r', 
                    figsize = 'auto',
                    save = None,
                    dot_size_norm = (5, 100),
                    dot_color_vmin = None,
                    dot_color_vmax = None,
                    show_plot = True,
                    comm_score_col = 'Commu_Score',
                    comm_score_cutoff = None,
                    cutoff_prop = None,
                    return_fig = False):
        """
        dot plot to show the summary of communication numbers between sender and receiver 

        Params
        -----
        pval_method
            should be one of ['zztest_pval', 'ttest_pval', 'ranksum_test_pval', 'permutation_test_pval', 'zztest_fdr', 'ttest_fdr', 'ranksum_test_fdr', 'permutation_test_fdr'], default is permutation_test_fdr
        pval_cutoff
            float, set to filter out non-significant communication events
        cmap
            color map to set dot color 
        figsize
            auto or a tuple of float such as (5.5, 4.2), defualt will be automatically estimate
        save
            str, the file name to save the figure
        dot_size_norm
            two values in a tuple, used to normalize the dot size, such as (10, 150)
        dot_color_vmin
            float, the value limits the color map in maximum
        dot_color_vmax
            float, the value limits the color map in minimum
        show_plot
            True or False, whether print the figure on the screen
        comm_score_col
            column name of communication score, can be Commu_Score
        comm_score_cutoff
            a float, set a cutoff so only communications with score greater than the cutoff will be focused
        cutoff_prop
            a float between 0 and 1, set a cutoff to filter out lowly abundant cell populations by the fraction of cells expressed sensor genes or metabolite, Note that this parameter will lost the function if cutoff_prop was set lower than the one user set at begaining of running mebocost.infer_commu or preparing mebocost object. This parameter were designed to further strengthen the filtering.
        return_fig:
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself
        """
        
        
#         if show_plot is None and self.show_plot is not None:
#             show_plot = self.show_plot

        comm_res = self.commu_res
        ## pdf
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None

        fig = CP._count_dot_plot_(commu_res=comm_res, pval_method = pval_method, pval_cutoff = pval_cutoff, 
                        cmap = cmap, figsize = figsize, pdf = Pdf, dot_size_norm = dot_size_norm, 
                        dot_color_vmin = dot_color_vmin, dot_color_vmax = dot_color_vmax, show_plot = show_plot,
                        comm_score_col = comm_score_col, comm_score_cutoff = comm_score_cutoff, cutoff_prop = cutoff_prop,
                        return_fig = return_fig)
        if save is not None and save is not False:
            Pdf.close()
        if return_fig:
            return(fig)

    def commu_network_plot(self,
                        sender_focus = [],
                        metabolite_focus = [],
                        sensor_focus = [],
                        receiver_focus = [],
                        remove_unrelevant = False,
                        and_or = 'and',
                        pval_method = 'permutation_test_fdr',
                        pval_cutoff = 0.05,
                        node_cmap = 'tab20',
                        figsize = 'auto',
                        line_cmap = 'RdBu_r',
                        line_color_vmin = None,
                        line_color_vmax = None,
                        linewidth_norm = (0.1, 1),
                        node_size_norm = (50, 300),
                        adjust_text_pos_node = True,
                        node_text_hidden = False,
                        node_text_font = 10,
                        save = None,
                        show_plot = True,
                        comm_score_col = 'Commu_Score',
                        comm_score_cutoff = None,
                        cutoff_prop = None,
                        text_outline = False,
                        return_fig = False):

        """
        Network plot to show the communications between cell groups

        Params
        ------
        sender_focus
            a list, set a list of sender cells to be focused, only plot related communications
        metabolite_focus
            a list, set a list of metabolites to be focused, only plot related communications
        sensor_focus
            a list, set a list of sensors to be focused, only plot related communications
        receiver_focus
            a list, set a list of receiver cells to be focused, only plot related communications
        remove_unrelevant
            True or False, set True to hide unrelated nodes
        and_or
            eithor 'and' or 'or', 'and' for finding communications that meet to all focus, 'or' for union
        pval_method
            should be one of ['zztest_pval', 'ttest_pval', 'ranksum_test_pval', 'permutation_test_pval', 'zztest_fdr', 'ttest_fdr', 'ranksum_test_fdr', 'permutation_test_fdr'], default is permutation_test_fdr
        pval_cutoff
            float, set to filter out non-significant communication events
        node_cmap
            node color map, used to indicate different cell groups, set one from https://matplotlib.org/stable/tutorials/colors/colormaps.html
        figsize
            auto or a tuple of float such as (5.5, 4.2), defualt will be automatically estimate
        line_cmap
            line color map, used to indicate number of communication events, set one from https://matplotlib.org/stable/tutorials/colors/colormaps.html
        line_color_vmin
            float, the value limits the line color map in minimum
        line_color_vmax
            float, the value limits the line color map in maximum
        linewidth_norm
            two values in a tuple, used to normalize the dot size, such as (0.1, 1)
        node_size_norm
            two values in a tuple, used to normalize the node size, such as (50, 300)
        adjust_text_pos_node 
            True or Flase, whether adjust the text position to avoid overlapping automatically
        node_text_font
            float, font size for node text annotaion
        save
            str, the file name to save the figure
        show_plot
            True or False, whether print the figure on the screen
        comm_score_col
            column name of communication score, can be Commu_Score
        comm_score_cutoff
            a float, set a cutoff so only communications with score greater than the cutoff will be focused
        cutoff_prop
            a float between 0 and 1, set a cutoff to filter out lowly abundant cell populations by the fraction of cells expressed sensor genes or metabolite, Note that this parameter will lost the function if cutoff_prop was set lower than the one user set at begaining of running mebocost.infer_commu or preparing mebocost object. This parameter were designed to further strengthen the filtering.
        return_fig:
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself
        """
        
        comm_res = self.commu_res
        ## pdf
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None

        fig = CP._commu_network_plot_(commu_res=comm_res, sender_focus = sender_focus, metabolite_focus = metabolite_focus, 
                            sensor_focus = sensor_focus, receiver_focus = receiver_focus, and_or = and_or, 
                            pval_method = pval_method, 
                            pval_cutoff = pval_cutoff, node_cmap = node_cmap, figsize = figsize, line_cmap = line_cmap, 
                            line_color_vmin = line_color_vmin, line_color_vmax = line_color_vmax,
                            linewidth_norm = linewidth_norm, node_text_hidden = node_text_hidden,
                            node_size_norm = node_size_norm, adjust_text_pos_node = adjust_text_pos_node, 
                            comm_score_col = comm_score_col, comm_score_cutoff = comm_score_cutoff, cutoff_prop = cutoff_prop,
                            node_text_font = node_text_font, pdf = Pdf, show_plot = show_plot, text_outline = text_outline,
                            return_fig = return_fig)
        if save is not None and save is not False:
            Pdf.close()
        
        if return_fig:
            return(fig)
            
    def violin_plot(self,
                    sensor_or_met,
                    cell_focus = [],
                    cmap = None,
                    vmin = None,
                    vmax = None,
                    figsize = 'auto',
                    cbar_title = '',
                    save = None,
                    show_plot = True,
                    return_fig = False):
        """
        Violin plot to show the distribution of sensor expression or metabolite level across cell groups

        Params
        -----
        sensor_or_met
            a list, provide a list of sensor gene name or metabolite name
        cell_focus
            a list, provide a list of cell type that you want to focus, otherwise keep empty
        cmap
            the color map used to draw the violin
        vmin
            float, maximum value for the color map
        vmin
            float, minimum value for the color map
        figsize
            auto or a tuple of float such as (5.5, 4.2), defualt will be automatically estimate
        title
            str, figure title on the top
        save
            str, the file name to save the figure
        show_plot
            True or False, whether print the figure on the screen
        comm_score_col
            column name of communication score, can be Commu_Score
        comm_score_cutoff
            a float, set a cutoff so only communications with score greater than the cutoff will be focused
        return_fig:
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself
        """
        ## cell group
        cell_ann = self.cell_ann.copy()
        if 'cell_group' not in cell_ann.columns.tolist():
            raise ValueError('ERROR: "cell_group" not in cell_ann column names!')
        ### extract expression for sensor
        sensors = []
        if self.exp_mat is not None and self.exp_mat_indexer is not None:
            sensor_loc = np.where(pd.Series(self.exp_mat_indexer).isin(sensor_or_met))
            #[i for i,j in enumerate(self.exp_mat_indexer.tolist()) if j in sensor_or_met]
            sensors = self.exp_mat_indexer[sensor_loc]
            #[j for i,j in enumerate(self.exp_mat_indexer.tolist()) if j in sensor_or_met]
            exp_dat = pd.DataFrame(self.exp_mat[sensor_loc].toarray(),
                                   index = sensors,
                                   columns = self.exp_mat_columns)
            if len(sensors) > 0:
                info('Find genes %s to plot violin'%(sensors))
                ## expression
                if save is not None and save is not False and isinstance(save, str):
                    save = save.replace('.pdf', '_sensor_exp.pdf')
                    Pdf = PdfPages(save)
                else:
                    Pdf = None
                if cmap is None:
                    ccmap = 'Reds'
                else:
                    ccmap = cmap

                if cbar_title == '':
                    sensor_cbar_title = 'Mean Expression'
                else:
                    sensor_cbar_title = cbar_title
                ## data mat for plot
                dat_mat = pd.merge(exp_dat.T, cell_ann[['cell_group']], left_index = True, right_index = True)
                fig = CP._violin_plot_(dat_mat=dat_mat, sensor_or_met=list(sensors),
                                       cell_focus = cell_focus, cmap = ccmap,
                                       vmin = vmin, vmax = vmax, figsize = figsize, 
                                       cbar_title = sensor_cbar_title, pdf = Pdf,
                                       show_plot = show_plot, return_fig = return_fig)

                if save is not None and save is not False:
                    Pdf.close()
                if return_fig:
                    return(fig)
            else:
                info('Warnings: no sensors to plot')
        else:
            info('Warnings: failed to load metabolite data matrix')
            
        ### extract metabolite level
        metabolites = list(set(sensor_or_met) - set(sensors))
        metabolites = list(set(metabolites) & set(self.met_ann['metabolite'].unique().tolist()))
        if metabolites:
            # to HMDBID
            met_name_to_id = {}
            for m, iD in self.met_ann[['metabolite', 'HMDB_ID']].values.tolist():
                met_name_to_id[m] = iD
            metaboliteIds = {x: met_name_to_id.get(x) for x in metabolites}
            ## metabolite matrix
            if self.met_mat is not None and self.met_mat_indexer is not None:
                met_loc = np.where(pd.Series(self.met_mat_indexer).isin(list(metaboliteIds.values())))[0]
                met_Ids = self.met_mat_indexer[met_loc]
                met_names = [list(metaboliteIds.keys())[list(metaboliteIds.values()).index(x)] for x in met_Ids]
                met_dat = pd.DataFrame(self.met_mat[met_loc].toarray(),
                                   index = met_names,
                                   columns = self.met_mat_columns)
                dat_mat = pd.merge(met_dat.T, cell_ann[['cell_group']], left_index = True, right_index = True)
                if len(met_names) > 0:
                    info("Find metabolites %s to plot violin"%metabolites)
                    ## expression
                    if save is not None and save is not False and isinstance(save, str):
                        save = save.replace('.pdf', '_metabolite.pdf')
                        Pdf = PdfPages(save)
                    else:
                        Pdf = None
                    if cmap is None:
                        ccmap = 'Purples'
                    else:
                        ccmap = cmap
                    if cbar_title == '':
                        met_cbar_title = 'Mean Abundance'
                    else:
                        met_cbar_title = cbar_title
                        
                    fig = CP._violin_plot_(dat_mat=dat_mat, sensor_or_met=list(metaboliteIds.keys()),
                                     cell_focus = cell_focus, cmap = ccmap,
                                    vmin = vmin, vmax = vmax, figsize = figsize,
                                    cbar_title = met_cbar_title, pdf = Pdf,
                                    show_plot = show_plot, return_fig = return_fig)

                    if save is not None and save is not False:
                        Pdf.close()
                    if return_fig:
                        return(fig)
                else:
                    info('Warnings: no metabolites to plot')
            else:
                info('Warnings: failed to load metabolite data matrix')
        else:
            info('Warnings: no metabolites to plot')
            
    def communication_in_notebook(self,
                                  pval_method = 'permutation_test_fdr',
                                  pval_cutoff = 0.05,
                                  comm_score_col = 'Commu_Score',
                                  comm_score_cutoff = None, 
                                  cutoff_prop = None
                                 ):

        # some handy functions to use along widgets
        from IPython.display import display, Markdown, clear_output, HTML
        import ipywidgets as widgets
        import functools

        outt = widgets.Output()

        df = self.commu_res.copy()
        
        if not comm_score_cutoff:
            comm_score_cutoff = 0
        if not cutoff_prop:
            cutoff_prop = 0
        ## basic filter
        df = df[(df[pval_method] <= pval_cutoff) & 
                (df[comm_score_col] >= comm_score_cutoff) &
                (df['metabolite_prop_in_sender'] >= cutoff_prop) &
                (df['sensor_prop_in_receiver'] >= cutoff_prop)
                ]
        
        senders = ['All']+sorted(list(df['Sender'].unique()))
        receivers = ['All']+sorted(list(df['Receiver'].unique()))
        metabolites = ['All']+sorted(list(df['Metabolite_Name'].unique()))
        transporters = ['All']+sorted(list(df['Sensor'].unique()))
        
        logic_butt = widgets.RadioButtons(
                            options=['and', 'or'],
                            description='Logic',
                            disabled=False
                        )

        sender_sel = widgets.SelectMultiple(description='Sender:',
                                            options=senders,
                                            layout=widgets.Layout(width='30%'))
        receiver_sel = widgets.SelectMultiple(description='Receiver:',
                                              options=receivers,
                                              layout=widgets.Layout(width='30%'))
        metabolite_sel = widgets.SelectMultiple(description='Metabolite:',
                                                options=metabolites,
                                                layout=widgets.Layout(width='30%'))
        sensor_sel = widgets.SelectMultiple(description='Sensor:',
                                                 options=transporters,
                                                layout=widgets.Layout(width='30%'))
        
        flux_butt = widgets.Button(description='Communication Flow (FlowPlot)',
                              layout=widgets.Layout(width='100%'))
        net_butt = widgets.Button(description='Communication Network (CirclePlot)',
                              layout=widgets.Layout(width='100%'))
        dotHeatmap_butt = widgets.Button(description='Communication Details (Dot-shaped Heatmap)',
                              layout=widgets.Layout(width='100%'))
        violin_butt = widgets.Button(description='ViolinPlot to show metabolite or sensor level in cell groups',
                              layout=widgets.Layout(width='100%'))

        def _flowplot_filter_(b):
            with outt:
                clear_output()
                print('+++++++++++++++++++++++++++ Running, Please Wait +++++++++++++++++++++++++++ ')
                print('[Selection]: Sender{}; Metabolite{}; Transporter{}; Receiver{}'.format(sender_sel.value,
                                                                                                  metabolite_sel.value,
                                                                                                  sensor_sel.value,
                                                                                                  receiver_sel.value))
                and_or = logic_butt.value
                
                self.FlowPlot(pval_method=pval_method,
                            pval_cutoff=pval_cutoff,
                            sender_focus = [x for x in sender_sel.value if x != 'All'],
                            metabolite_focus = [x for x in metabolite_sel.value if x != 'All'],
                            sensor_focus = [x for x in sensor_sel.value if x != 'All'],
                            receiver_focus = [x for x in receiver_sel.value if x != 'All'],
                            remove_unrelevant = True,
                            and_or = and_or,
                            node_label_size = 8,
                            node_alpha = .8,
                            figsize = 'auto',
                            node_cmap = 'Set1',
                            line_cmap = 'bwr',
                            line_vmin = None,
                            line_vmax = None,
                            node_size_norm = (10, 150),
                            linewidth_norm = (0.5, 5),
                            save=None, 
                            show_plot = False,
                            comm_score_col = comm_score_col,
                            comm_score_cutoff = comm_score_cutoff,
                            cutoff_prop = cutoff_prop,
                            text_outline = False,
                            return_fig = False)
                
                
        def _networkplot_filter_(b):
            with outt:
                clear_output()
                print('+++++++++++++++++++++++++++ Running, Please Wait +++++++++++++++++++++++++++ ')
                print('[Selection]: Sender{}; Metabolite{}; Transporter{}; Receiver{}'.format(sender_sel.value,
                                                                                                  metabolite_sel.value,
                                                                                                  sensor_sel.value,
                                                                                                  receiver_sel.value))
                and_or = logic_butt.value
                self.commu_network_plot(
                                sender_focus = [x for x in sender_sel.value if x != 'All'],
                                metabolite_focus = [x for x in metabolite_sel.value if x != 'All'],
                                sensor_focus = [x for x in sensor_sel.value if x != 'All'],
                                receiver_focus = [x for x in receiver_sel.value if x != 'All'],
                                remove_unrelevant = False,
                                and_or = and_or,
                                pval_method = pval_method,
                                pval_cutoff = pval_cutoff,
                                node_cmap = 'tab20',
                                figsize = 'auto',
                                line_cmap = 'RdBu_r',
                                line_color_vmin = None,
                                line_color_vmax = None,
                                linewidth_norm = (0.1, 1),
                                node_size_norm = (50, 300),
                                adjust_text_pos_node = False,
                                node_text_font = 10,
                                save = None,
                                show_plot = True,
                                comm_score_col = comm_score_col,
                                comm_score_cutoff = comm_score_cutoff,
                                cutoff_prop = cutoff_prop,
                                text_outline = False
                                )
        def _dotHeatmapPlot_(b):
            with outt:
                clear_output()
                print('+++++++++++++++++++++++++++ Running, Please Wait +++++++++++++++++++++++++++ ')
                print('[Selection]: Sender{}; Metabolite{}; Transporter{}; Receiver{}'.format(sender_sel.value,
                                                                                                  metabolite_sel.value,
                                                                                                  sensor_sel.value,
                                                                                                  receiver_sel.value))
                and_or = logic_butt.value
                self.commu_dotmap(
                            sender_focus = [x for x in sender_sel.value if x != 'All'],
                            metabolite_focus = [x for x in metabolite_sel.value if x != 'All'],
                            sensor_focus = [x for x in sensor_sel.value if x != 'All'],
                            receiver_focus = [x for x in receiver_sel.value if x != 'All'],
                            and_or = and_or,
                            pval_method=pval_method,
                            pval_cutoff=pval_cutoff, 
                            figsize = 'auto',
                            cmap = 'bwr',
                            node_size_norm = (10, 150),
                            save = None, 
                            show_plot = True,
                            comm_score_col = comm_score_col,
                            comm_score_cutoff = comm_score_cutoff,
                            cutoff_prop = cutoff_prop
                )

        def _violinPlot_(b):
            with outt:
                clear_output()
                print('+++++++++++++++++++++++++++ Running, Please Wait +++++++++++++++++++++++++++ ')
                print('[Selection]: Sender{}; Metabolite{}; Transporter{}; Receiver{}'.format(sender_sel.value,
                                                                                                  metabolite_sel.value,
                                                                                                  sensor_sel.value,
                                                                                                  receiver_sel.value))
                
                self.violin_plot(
                                sensor_or_met = [x for x in metabolite_sel.value + sensor_sel.value if x != 'All'],
                                cell_focus = [x for x in sender_sel.value + receiver_sel.value if x != 'All'],
                                cmap = None,
                                vmin = None,
                                vmax = None,
                                figsize = 'auto',
                                cbar_title = '',
                                save = None,
                                show_plot = True)
                
                
        flux_butt.on_click(_flowplot_filter_)
        net_butt.on_click(_networkplot_filter_)
        dotHeatmap_butt.on_click(_dotHeatmapPlot_)
        violin_butt.on_click(_violinPlot_)


        h1 = widgets.HBox([sender_sel, metabolite_sel, sensor_sel, receiver_sel])
        h2 = widgets.VBox([flux_butt, net_butt, dotHeatmap_butt, violin_butt])

        mk = Markdown("""<b>Select and Click button to visulize</b>""")
        display(mk, widgets.VBox([logic_butt, h1, h2, outt]))

## ===================== pathway plot functions =============
        
    def pathway_scatter(self, 
                a_pair, 
                pval_cutoff=0.05, 
                ES_cutoff=0,
                cmap = 'cool',
                vmax = None,
                vmin = None,
                figsize = 'auto',
                title = '',
                maxSize = 500,
                minSize = 15,
                save = None,
                show_plot = True,
                return_fig = False):
        """
        Plot the associated pathway in scatter plot
        
        Params
        -----
        a_pair
            str, the format should be either sensor ~ receiver or sender ~ receiver
        pval_cutoff
            float, cutoff of pval (padj) to focus on significantly associated pathways
        ES_cutoff
            float, cutoff of Fold Enrichment Score for assiciated pathways, positive value for positively associated pathways
        cmap
            python color map for showing significance
        vmax
            the maximum limits for colormap
        vmin
            the minimum limits for colormap
        figsize
            a tuple inclues two values, represents width and height, set "auto" to automatically estimate
        title:
            str, figure title
        maxSize
            the term size in maximum
        minSize
            the term size in mimum
        save
            str, the path of where the figure save to
        show_plot
            True or False, to display the figure on the screen or not
        return_fig:
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself
        """
        
        cellpairs = list(self.enrich_result['cellpair_res'].keys())
        sensor_receivers = list(self.enrich_result['sensor_res'].keys())

        if a_pair in cellpairs:
            res_dict = self.enrich_result['cellpair_res']
        if a_pair in sensor_receivers:
            res_dict = self.enrich_result['sensor_res']
        if a_pair not in cellpairs and a_pair not in sensor_receivers:
            raise KeyError('ERROR to read given a_pair!')
        
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None

        fig = PP._scatter_(res_dict = res_dict, a_pair = a_pair, pval_cutoff=pval_cutoff, 
                    ES_cutoff=ES_cutoff, cmap = cmap, vmax = vmax, vmin = vmin,
                    figsize = figsize, title = title, maxSize = maxSize, minSize = minSize,
                    pdf = Pdf, show_plot = show_plot, return_fig = return_fig)
        
        if save is not None and save is not False:
            Pdf.close()
        if return_fig:
            return(fig)
        
    def pathway_stacked_bar(self,
                pair1,
                pair2,
                pval_cutoff=0.05, 
                ES_cutoff=0,
                cmap = 'spring_r',
                vmax = None,
                vmin = None,
                figsize = 'auto',
                title = '',
                maxSize = 500,
                minSize = 15,
                colors = ['#CC6677', '#1E90FF'],
                save = None,
                show_plot = True,
                return_fig = False):
    
        """
        compare pathways for two communication events
        
        Params
        ------
        pair1 and pair2
            str, the format of pair1 and pairs should be both sensor ~ receiver, or both sender ~ receiver
        pval_cutoff
            float, cutoff of pval (padj) to focus on significantly associated pathways
        ES_cutoff
            float, cutoff of Fold Enrichment Score for assiciated pathways, positive value for positively associated pathways
        cmap
            python color map for showing significance
        vmax
            the maximum limits for colormap
        vmin
            the minimum limits for colormap
        figsize
            a tuple inclues two values, represents width and height, set "auto" to automatically estimate
        title:
            str, figure title
        maxSize
            the term size in maximum
        minSize
            the term size in mimum
        save
            str, the path of where the figure save to
        show_plot
            True or False, to display the figure on the screen or not
        return_fig
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself
        """
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None
        
        cellpairs = list(self.enrich_result['cellpair_res'].keys())
        sensor_receivers = list(self.enrich_result['sensor_res'].keys())

        cellpair = np.all([pair1 in cellpairs, pair2 in cellpairs])
        sensor_receiver = np.all([pair1 in sensor_receivers, pair2 in sensor_receivers])

        if cellpair:
            res_dict = self.enrich_result['cellpair_res']
        elif sensor_receiver:
            res_dict = self.enrich_result['sensor_res']
        else:
            raise KeyError('pair1 and pair2 should be both sender-receiver or sensor-receiver, please check!')
        
            
        fig = PP._stacked_bar_(res_dict = res_dict, pair1 = pair1, pair2 = pair2, 
                         pval_cutoff = pval_cutoff, ES_cutoff = ES_cutoff,
                         cmap = cmap, vmax = vmax, vmin = vmin, figsize = figsize,
                         title = title, maxSize = maxSize, minSize = minSize, colors = colors,
                         pdf = Pdf, show_plot = show_plot, return_fig = return_fig)
                
        if save is not None and save is not False:
            Pdf.close()
        if return_fig:
            return(fig)
            
    
    def pathway_multi_dot(self,
                    pairs, 
                    pval_cutoff=0.05, 
                    ES_cutoff=0,
                    cmap = 'Spectral_r',
                    vmax = None,
                    vmin = None,
                    node_size_norm = (20, 100),
                    figsize = 'auto',
                    title = '',
                    maxSize = 500,
                    minSize = 15,
                    save = None,
                    show_plot = True,
                    swap_axis = False,
                    return_fig = False):
        """
        draw dot map to show the associated pathways in multiple comparisons
        
        Params
        -----
        pairs
            a list, elements should be all in sender-receiver or all in sensor-receiver
        pval_cutoff
            float, cutoff of pval (padj) to focus on significantly associated pathways
        ES_cutoff
            float, cutoff of Normalized Enrichment Score for assiciated pathways, positive value for positively associated pathways, negative value for negatively associated
        cmap
            python color map for showing significance
        vmax
            the maximum limits for colormap
        vmin
            the minimum limits for colormap
        figsize
            a tuple inclues two values, represents width and height, set "auto" to automatically estimate
        title:
            str, figure title
        maxSize
            the term size in maximum
        minSize
            the term size in mimum
        save
            str, the path of where the figure save to
        show_plot
            True or False, to display the figure on the screen or not
        swap_axis
           True or False, set True to flip x and y axis
        return_fig
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself
        """
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None
        
        cellpairs = list(self.enrich_result['cellpair_res'].keys())
        sensor_receivers = list(self.enrich_result['sensor_res'].keys())
        
        if np.all([x in cellpairs for x in pairs]):
            res_dict = self.enrich_result['cellpair_res']
        elif np.all([x in sensor_receivers for x in pairs]):
            res_dict = self.enrich_result['sensor_res']
        else:
            raise KeyError('pairs should be a list of elements all in sender-receiver or all in sensor-receiver, please check!')

        fig = PP._multi_dot_(res_dict = res_dict, pairs = pairs, pval_cutoff=pval_cutoff, ES_cutoff=ES_cutoff,
                    cmap = cmap, vmax = vmax, vmin = vmin, node_size_norm = node_size_norm,
                    figsize = figsize, title = title, maxSize = maxSize, minSize = minSize, pdf = Pdf,
                    show_plot = show_plot, swap_axis = swap_axis, return_fig = return_fig)
        
        if save is not None and save is not False:
            Pdf.close()    
        if return_fig:
            return(fig)
        
    def pathway_ES_plot(self, 
                  a_pair, 
                  pathway,
                  figsize = (8, 3.5),
                  dot_color = '#1874CD',
                  curve_color = 'black',
                  title='',
                  save = None,
                  show_plot = True,
                  return_fig = False,
                  return_data = False
                 ):
        """
        this function to draw enrichment plot for certain pathway in certain cell
        
        Params
        ------
        a_pair
            str, the format should be either sensor ~ receiver or sender ~ receiver
        pathway
            str, pathway name
        figsize
            a tuple with two values, use to set the figure size, default is (10, 5)
        dot_color
            color name, use to set the pathway genes in the scatter plot
        line_color
            color name, use to set the pathway genes in enrichement plot
        ht_cmap
            color map, use to set the color scale for background
        title
            str, figure title
        save
            path, save the figure to
        show_plot
            True or False, display the figure on the screen or not
        return_fig
            True or False, set True to return the figure object, this can be useful if you want to manipulate figure by yourself        
        """
        if save is not None and save is not False and isinstance(save, str):
            Pdf = PdfPages(save)
        else:
            Pdf = None
        
        ## geneset
        all_pathway = self.enrich_result['gmt']
        
        ## remove kegg id
        if not re.search('hsa[0-9]* |mmu[0-9]* ', pathway):
            all_pathway_update = {k.replace(re.search('hsa[0-9]* |mmu[0-9]* ', k).group(), '') : all_pathway[k] for k in all_pathway.keys()}
            name_match = {k.replace(re.search('hsa[0-9]* |mmu[0-9]* ', k).group(), '') : k for k in all_pathway.keys()}
        if pathway not in list(all_pathway_update.keys()):
            raise KeyError('ERROR: cannot find the pathway in database, please check the name!')
        geneSet = all_pathway_update[pathway]
        pathway = name_match[pathway]
        ##
        cellpairs = list(self.enrich_result['cellpair_res'].keys())
        sensor_receivers = list(self.enrich_result['sensor_res'].keys())

        ## sensor in receiver
        if a_pair in sensor_receivers:
            s, r = a_pair.split(' ~ ')
            s_loc = np.where(self.enrich_result['weighted_exp'][s]['index'] == r)
            weightList = pd.Series(self.enrich_result['weighted_exp'][s]['weight_exp'][s_loc].toarray()[0],
                                   index = self.enrich_result['weighted_exp'][s]['columns'])
            expList = self.enrich_result['avg_exp_norm'][r]
            gxtList = self.gene_network[s]
            mHG_obj = self.enrich_result['sensor_res'][a_pair]['mHG_obj'][pathway]
            
            fig = PP._sensor_ES_plot_(geneSet=geneSet, 
                          mHG_obj=mHG_obj,
                          expList=expList,
                          gxtList=gxtList,
                          sensor = s,
                          receiver = r,
                          figsize = figsize,
                          dot_color = dot_color,
                          curve_color = curve_color,
                          title = title,
                          pdf = Pdf,
                          show_plot = show_plot,
                          return_fig = return_fig
                     )

        ## sender to receiver
        if a_pair in cellpairs:
            sender, receiver = a_pair.split(' ~ ')
            weightList = self.enrich_result['weighted_exp_agg'][a_pair]
            expList = self.enrich_result['avg_exp_norm'][receiver]
            mHG_obj = self.enrich_result['cellpair_res'][a_pair]['mHG_obj'][pathway]

            fig = PP._cellpair_ES_plot_(geneSet=geneSet, 
                      mHG_obj=mHG_obj,
                      expList=expList,
                      gxtList=weightList,
                      sender = sender,
                      receiver = receiver,
                      figsize = figsize,
                      dot_color = dot_color,
                      curve_color = curve_color,
                      title = title,
                      pdf = Pdf,
                      show_plot = show_plot,
                      return_fig = return_fig
                 )
            
                
        if a_pair not in cellpairs and a_pair not in sensor_receivers:
            raise KeyError('ERROR to read given a_pair!')
        
        if save is not None and save is not False:
            Pdf.close()    
        
        if return_data and return_fig:
            df = pd.concat([expList, gxtList], axis = 1).reset_index().dropna()
            df.columns = ['gene', 'scaled_expression', 'gene_correlation']
            df = df.loc[df['gene'].isin(geneSet)]
            return(fig, df)
        elif return_data and not return_fig:
            df = pd.concat([expList, weightList], axis = 1).reset_index().dropna()
            df.columns = ['gene', 'scaled_expression', 'gene_weight']
            df = df.loc[df['gene'].isin(geneSet)]
            return(df)
        elif return_fig and not return_data:
            return(fig)
        else:
            pass
     
    def pathway_in_notebook(self):

        # some handy functions to use along widgets
        from IPython.display import display, Markdown, clear_output, HTML
        # widget packages
        import ipywidgets as widgets
        import functools

        
        ## sensor ~ cell
        sender_receiver = list(self.enrich_result['sensor_res'].keys())
        sender_receiver = sorted(list(set(sender_receiver)))

        ## cell -> cell
        cellpair = list(self.enrich_result['cellpair_res'].keys())
        cellpair = sorted(list(set(cellpair)))
        
        all_pathway = self.enrich_result['gmt']
        all_pathway = {k.replace(re.search('hsa[0-9]* |mmu[0-9]* ', k).group(), '') : all_pathway[k] for k in all_pathway.keys()}
        all_pathway_names = list(sorted(list(all_pathway.keys())))

        # creating menu with them 
        pathway_vars = widgets.Select(
                        description = 'Pathway:',
                        options=all_pathway_names,
                        layout={'width': '90%'}, # If the items' names are long
                        disabled=False,
                    )
        
        sender_receiver_sel = widgets.SelectMultiple(options=sender_receiver,
                              layout=widgets.Layout(width='90%'))
        cellpair_sel = widgets.SelectMultiple(options=cellpair,
                              layout=widgets.Layout(width='90%'))

        # button, output, function and linkage
        senser_receiver_butt = widgets.Button(description='CLICK to Print Pathway for Sensor~Receiver Cell',
                              layout=widgets.Layout(width='90%'))
        cellpair_butt = widgets.Button(description='CLICK to Print Pathway for Sender Cell -> Receiver Cell',
                              layout=widgets.Layout(width='90%'))
        senser_receiver_pathway_butt = widgets.Button(description='CLICK to show enrichment curve for Sensor~Receiver Cell',
                              layout=widgets.Layout(width='90%'))
        cellpair_pathway_butt = widgets.Button(description='CLICK to show enrichment curve for Sender Cell -> Receiver Cell',
                              layout=widgets.Layout(width='90%'))
        
        outt = widgets.Output()

        def _one_single_clicked(b, Type):
            with outt:
                clear_output()
                print('WARINING: Running, Please do not reflesh until you see the figure!')
                pairs = sender_receiver_sel.value if Type == 'sensor' else cellpair_sel.value
                if len(pairs) == 1:
                    self.pathway_scatter( 
                            a_pair = pairs[0], 
                            pval_cutoff=0.05, 
                            ES_cutoff=0,
                            cmap = 'cool',
                            vmax = None,
                            vmin = None,
                            figsize = 'auto',
                            title = '',
                            maxSize = 500,
                            minSize = 15,
                            save = None,
                            show_plot = True)

                elif len(pairs) == 2:
                    self.pathway_stacked_bar(
                                pair1 = pairs[0],
                                pair2 = pairs[1],
                                pval_cutoff=0.05, 
                                ES_cutoff=0,
                                cmap = 'spring_r',
                                vmax = None,
                                vmin = None,
                                figsize = 'auto',
                                title = '',
                                maxSize = 500,
                                minSize = 15,
                                colors = ['#CC6677', '#1E90FF'],
                                save = None,
                                show_plot = True)
                else:
                    self.pathway_multi_dot(
                                pairs = pairs, 
                                pval_cutoff=0.05, 
                                ES_cutoff=0,
                                cmap = 'Spectral_r',
                                vmax = None,
                                vmin = None,
                                node_size_norm = (20, 100),
                                figsize = 'auto',
                                title = '',
                                maxSize = 500,
                                minSize = 15,
                                save = None,
                                show_plot = True)


        def _enrich_clicked(b, Type):
            with outt:
                clear_output()
                print('WARINING: Running, Please do not reflesh until you see the figure!')

                pairs = sender_receiver_sel.value if Type == 'sensor' else cellpair_sel.value
                a_pair = pairs[0] # only one used
                print(a_pair)
                ## pathway value
                pathway = pathway_vars.value
                print(pathway)
                ## extract ES and padj
                ##
                cellpairs = list(self.enrich_result['cellpair_res'].keys())
                sensor_receivers = list(self.enrich_result['sensor_res'].keys())

                ## plot
                self.pathway_ES_plot(a_pair = a_pair, 
                              pathway = pathway,
                              figsize = (8, 3.5),
                              dot_color = '#1874CD',
                              curve_color = 'black',
                              title='',
                              save = None,
                              show_plot = True,
                              return_fig = False,
                              return_data = False
                             )
               
        
        senser_receiver_butt.on_click(functools.partial(_one_single_clicked, Type = 'sensor'))
        # display for sensor in receiver
        box1 = widgets.VBox([sender_receiver_sel, senser_receiver_butt],
                            layout=widgets.Layout(width='50%'))

        cellpair_butt.on_click(functools.partial(_one_single_clicked, Type = 'cellpair'))
        # display for sender to receiver
        box2 = widgets.VBox([cellpair_sel, cellpair_butt],
                            layout=widgets.Layout(width='50%'))
        
        box3 = widgets.VBox([pathway_vars,
                             widgets.HBox([senser_receiver_pathway_butt,
                                           cellpair_pathway_butt])])
        senser_receiver_pathway_butt.on_click(functools.partial(_enrich_clicked, Type = 'sensor'))
        cellpair_pathway_butt.on_click(functools.partial(_enrich_clicked, Type = 'cellpair'))

        mk = Markdown("""<b>Select one or multiple to visulize</b>""")
        display(mk, widgets.VBox([widgets.HBox([box1, box2]), box3, outt]))

            
            
            
            
            
            
            
            
            
            
