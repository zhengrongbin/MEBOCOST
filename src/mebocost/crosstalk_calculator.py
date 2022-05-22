#!/usr/bin/env python

# ================================
# @auther: Rongbin Zheng
# @email: Rongbin.Zheng@childrens.harvard.edu
# @date: May 2022
# ================================

import os,sys
import time
import pickle as pk
from datetime import datetime
import numpy as np
import pandas as pd
import traceback
from operator import itemgetter
import statsmodels
import statsmodels.api as sm
import scipy
from scipy import sparse
import collections
import multiprocessing

"""
main functions of metabolic cross talk  

"""


def info(string):
    """
    print information
    """
    today = datetime.today().strftime("%B %d, %Y")
    now = datetime.now().strftime("%H:%M:%S")
    current_time = today + ' ' + now
    print("[{}]: {}".format(current_time, string))


class InferComm:
    """
    class for infering communication
    input at least includes:
    expression data of metabolite sensor
    estimated metabolite level in cell-wise
    metabolite sensor inofrmation
    cell annotation used for grouping cells 

    Params
    -----
    exp_mat
        data frame, expression matrix, cells in row names, genes in column names
    met_mat
        data frame, metabolite abundance matrix, cells in row names, metabolite in column names
    cell_ann
        data frame, cell annotation information, cells in row names, columns are needed annotation such as cell_group, cell_type, cluster,
        by default, we automatically find "cell_type" and "cluster" columns to group cells and take the average sensor expression and metabolite abundance for them
    met_sensor
        data frame, metabolite sensor information, each row is a pair of metabolite and sensor, must include columns named HMDB_ID and Gene_name,
        we will model the communication for each pair of metabolite and sensor
    group_col
        a list, specify the column names in 'cell_ann' for grouping cells, by default 'cell_type' and 'cluster' will be used
    sensor_type
        a list, provide a list of sensor type that will be used in the communication modeling, must be one or more from ['receptor', 'transporter', 'interacting'], default is all the three

    """
    def __init__(self, 
                 exp_mat, 
                 exp_mat_indexer, 
                 exp_mat_columns,
                 met_mat,
                 met_mat_indexer, 
                 met_mat_columns,
                 cell_ann,
                 met_sensor,
                 avg_exp = None,
                 avg_exp_indexer = None,
                 avg_exp_columns = None,
                 avg_met = None,
                 avg_met_indexer = None,
                 avg_met_columns = None,
                 sensor_type=['Receptor', 'Transporter', 'Nuclear Receptor'], 
                 thread=None
                ):
        self.exp_mat = exp_mat
        self.exp_mat_indexer = exp_mat_indexer
        self.exp_mat_columns = exp_mat_columns

        self.met_mat = met_mat
        self.met_mat_indexer = met_mat_indexer
        self.met_mat_columns = met_mat_columns

        self.cell_ann = cell_ann
        
        self.avg_exp = avg_exp
        self.avg_exp_indexer = avg_exp_indexer
        self.avg_exp_columns = avg_exp_columns
        
        self.avg_met = avg_met
        self.avg_met_indexer = avg_met_indexer
        self.avg_met_columns = avg_met_columns
        
        ## focus on the given sensor types
        met_sensor = met_sensor[met_sensor['Annotation'].str.upper().isin([x.upper() for x in sensor_type])]
        ## only focus on met and gene in the data
        met_sensor = met_sensor[met_sensor['Gene_name'].isin(self.exp_mat_indexer) &
                                met_sensor['HMDB_ID'].isin(self.met_mat_indexer)]
        info('Sensor type used %s'%met_sensor['Annotation'].unique().tolist())
        self.met_sensor = met_sensor
        self.thread = thread
    

    def _get_shuffled_avg_exp_(self, i):
        """
        get averaged expression of met and gene in each permutation
        given i is a index of permutation
        """
        indexer = self.shuffle_index[i]
        cell_label_new = np.array(self.cell_ann['cell_group'].tolist())[indexer] ## shuffling the cell labels

        e_avg = np.empty(shape = (self.exp_mat.shape[0],0)) ## save exp data
        m_avg = np.empty(shape = (self.met_mat.shape[0],0)) ## save met data

        for x in self.group_names:
            cell_indexer = np.where(cell_label_new == x)[0]
            e_avg = np.concatenate((e_avg, self.exp_mat[:,cell_indexer].mean(axis = 1)), axis = 1)
            m_avg = np.concatenate((m_avg, self.met_mat[:,cell_indexer].mean(axis = 1)), axis = 1)
        
        e_avg = sparse.csr_matrix(e_avg)
        m_avg = sparse.csr_matrix(m_avg)
        return(i, e_avg, m_avg)
        
    
    def _shuffling_(self, group_names, n_shuffle = 1000, seed = 12345):
        """
        permutation to generate random backgroud
        """
        ## intersect cells in exp_mat, met_mat, and cell_ann
        common_cell = list(set(self.exp_mat_columns.tolist()) & 
                          set(self.met_mat_columns.tolist()) &
                          set(self.cell_ann.index.tolist()))
        ## cell ann
        self.cell_ann = self.cell_ann.loc[common_cell,:]
        ## exp
        cell_index = np.where(pd.Series(self.exp_mat_columns).isin(common_cell))[0]
        self.exp_mat = self.exp_mat[:,cell_index]
        self.exp_mat_columns = self.exp_mat_columns[cell_index]
        ## met
        cell_index = np.where(pd.Series(self.met_mat_columns).isin(common_cell))[0]
        self.met_mat = self.met_mat[:,cell_index]
        self.met_mat_columns = self.met_mat_columns[cell_index]
        
        ## get index of cells in shuffling
        np.random.seed(seed) ## set seed for reproducible shuffling
        self.shuffle_index = [np.random.choice(range(len(common_cell)),
                                               len(common_cell),
                                               replace = False) for i in range(n_shuffle)]
        
        ## cell group used to take mean
        cell_group = self.cell_ann.reindex(index = common_cell)[['cell_group']]
        ## exp mat for each shuffle
        info('take exp and met avg for shuffling')
        sensor_gene_loc = np.where(pd.Series(self.exp_mat_indexer).isin(self.met_sensor['Gene_name']))[0]
        sensor_gene = self.exp_mat_indexer[sensor_gene_loc]

        self.exp_mat = self.exp_mat[sensor_gene_loc,:]
        self.exp_mat_indexer = self.exp_mat_indexer[sensor_gene_loc]
        
        self.perm_exp_avg = collections.defaultdict()
        self.perm_met_avg = collections.defaultdict()
        
        pool = multiprocessing.Pool(self.thread)
        perm_col = pool.map(self._get_shuffled_avg_exp_, range(len(self.shuffle_index)))
        pool.close()
        ## into dict
        for i, e, m in perm_col:
            ## i is index of permute, e is exp shuffling, m is met shuffling
            self.perm_exp_avg[i] = e
            self.perm_met_avg[i] = m

        self.permute_exp_indexer = np.array(sensor_gene)
        self.permute_exp_columns = self.group_names
        
        self.permute_met_indexer = self.met_mat_indexer
        self.permute_met_columns = self.group_names
            
  
    def _commu_score_for_one_(self, 
                              avg_exp_sensor,
                              avg_met_met, 
                              perm_exp_avg_sensor,
                              perm_met_avg_met,
                              norm = None,
                              method = 'product'
                             ):
        """
        given a pair of sensor and metabolite
        calculate the communication score for eal data using ave_exp and avg_met object
        and calculate the backgroud of communication score using perm_exp_avg and perm_met_avg
        """
        ## real data
        sx = avg_exp_sensor.copy()
        mx = avg_met_met.copy()

        ## caculating communication score by products in pairwise cell group
        info('Normalizing Cluster Mean and Calculating Communication Score')
        norm_func = lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)) if np.max(x) != np.min(x) else np.array([0]*len(x))
        if norm == 'min-max':
            sx = norm_func(sx)
            mx = norm_func(mx)
            
        prod_mat = collections.defaultdict()
        if method == 'product':
            for i in range(len(self.group_names)):
                v = mx[i]
                c = self.group_names[i]
                ## met in one cell groups multiply all cell groups for sensor, so do pairwise
                prod_mat[c] = np.multiply(v, sx)
        elif method == 'average':
            for i in range(len(self.group_names)):
                v = mx[i]
                c = self.group_names[i]
                prod_mat[c] = np.array([np.mean([x, v]) for x in sx.tolist()])
        elif method == 'product_sqrt':
            for i in range(len(self.group_names)):
                v = mx[i]
                c = self.group_names[i]
                prod_mat[c] = np.nan_to_num(np.sqrt(np.multiply(v, sx)))
        else:
            raise KeyError('scoring problem!')

        prod_mat = pd.DataFrame.from_dict(prod_mat, orient = 'columns') ## column = sender
        prod_mat.index = self.group_names ## receiver
        ## reshape matrix
        prod_mat_values = prod_mat.unstack().reset_index()
#         prod_mat_values = prod_mat.where(np.triu(np.ones(prod_mat.shape)).astype(np.bool_)).unstack().reset_index().dropna()
        prod_mat_values.columns = ['Sender', 'Receiver', 'Commu_Score']

        ## background using shuffling data
        background = []
        if norm == 'min-max':
            ## norm by max-min normalization
            perm_exp_avg_sensor = {i: norm_func(perm_exp_avg_sensor[i]) for i in perm_exp_avg_sensor}
            perm_met_avg_met = {i: norm_func(perm_met_avg_met[i]) for i in perm_met_avg_met}
            
        if method == 'product':
            for i in range(len(perm_exp_avg_sensor)):
                e, m = perm_exp_avg_sensor[i], perm_met_avg_met[i]
                background.extend([[i, self.group_names[x]]+np.multiply(m[x], e).tolist() for x in range(len(self.group_names))])
        elif method == 'average':
            for i in range(len(perm_exp_avg_sensor)):
                e, m = perm_exp_avg_sensor[i], perm_met_avg_met[i]
                background.extend([[i, self.group_names[x]]+[np.mean([m[x], e]) for i in e.tolist()] for x in range(len(self.group_names))])
        elif method == 'product_sqrt':
            for i in range(len(perm_exp_avg_sensor)):
                e, m = perm_exp_avg_sensor[i], perm_met_avg_met[i]
                background.extend([[i, self.group_names[x]]+np.nan_to_num(np.sqrt(np.multiply(m[x], e))).tolist() for x in range(len(self.group_names))])
        else:
            raise KeyError('scoring problem!')

        background_df = pd.DataFrame(background, 
                          columns = ['N_permut', 'Sender_']+self.group_names)

        return prod_mat_values, background_df

    
    def cummin(self, x):
        """A python implementation of the cummin function in R"""
        for i in range(1, len(x)):
            if x[i-1] < x[i]:
                x[i] = x[i-1]
        return x


    def bh_fdr(self, pval):
        """A python implementation of the Benjamani-Hochberg FDR method.
        This code should always give precisely the same answer as using
        p.adjust(pval, method="BH") in R.
        Parameters
        ----------
        pval : list or array
            list/array of p-values
        Returns
        -------
        pval_adj : np.array
            adjusted p-values according the benjamani-hochberg method
        """
        pval_array = np.array(pval)
        sorted_order = np.argsort(pval_array)
        original_order = np.argsort(sorted_order)
        pval_array = pval_array[sorted_order]

        # calculate the needed alpha
        n = float(len(pval))
        pval_adj = np.zeros(int(n))
        i = np.arange(1, int(n)+1, dtype=float)[::-1]  # largest to smallest
        pval_adj = np.minimum(1, self.cummin(n/i * pval_array[::-1]))[::-1]
        return pval_adj[original_order]

    def _fdr_(self, 
              pvalue_res,
              testing_method = ['ztest', 'ttest', 'ranksum_test', 'permutation_test']):
        """
        calculate FDR correction
        """
        pvalue_res = pvalue_res.fillna(1) ## if NA
        for method in testing_method:
            fdr = self.bh_fdr(pvalue_res[method+'_pval'])
            pvalue_res[method+'_fdr'] = fdr
            del fdr
        return(pvalue_res)


    def _testing_(self, real_score, bg_values, method = 'ztest'):
        """
        run hypothesis testing
        real_score: the observed value
        bg_value: a array-like or list of background value
        """
        stat, pval = 0, 1
        if method == 'ztest':
            stat, pval = statsmodels.api.stats.ztest(bg_values,
                                    value = real_score, 
                                    alternative = "smaller")
        elif method == 'ttest':
            stat, pval, df = statsmodels.api.stats.ttest_ind(x1 = bg_values,
                                    x2 = [real_score], 
                                    alternative = "smaller")
        elif method == 'ranksum_test':
            stat, pval = scipy.stats.ranksums(x = bg_values, y = [real_score], alternative='less')
        # elif method == 'sign_test':
        #     stat, pval = statsmodels.stats.descriptivestats.sign_test(bg_values, real_score)
        elif method == 'permutation_test':
            bg_larger = 0
            for i in bg_values:
                if i > real_score:
                    bg_larger += 1
            stat, pval = bg_larger, bg_larger / len(bg_values)
        else:
            info('+++ unknown testing')
        return(stat, pval)

    def _signif_test_(self, commu_mat, background_df):
        """
        different statistical testing was used to evaluate the significance of observed communication score in background
        ---------
        """
        info('Calculating P-value')
        testing_method = ['ztest', 'ttest', 'ranksum_test', 'permutation_test']
        pvalue_res = []
        for n, iteract in commu_mat.sort_values('Commu_Score', ascending = False).iterrows():
            sender, receiver, real_score = iteract.tolist()
            bg_values = background_df[background_df['Sender_']==sender][receiver]
            ## z test to be used to compare population means to a sample's
            ## null hypothesis: mean of background is far more samller than observation
            testing_res = iteract.copy()
            testing_res['bg_mean'] = bg_values.mean()
            testing_res['bg_std'] = bg_values.std()
#             if testing_res['bg_mean'] != 0:
#                 testing_res['Norm_Commu_Score'] = testing_res['Commu_Score'] / testing_res['bg_mean']
#             else:
#                 testing_res['Norm_Commu_Score'] = testing_res['Commu_Score']
            ## hypothesis testing using more methods
            for method in testing_method:
                if testing_res['bg_std'] == 0 and testing_res['bg_mean'] == 0:
                    ## sometimes, the gene is not expressed at all cells, so no need to do testing
                    stat, pval = 0, 1
                else:
                    stat, pval = self._testing_(real_score=real_score, bg_values=bg_values, method = method)
                testing_res[method+'_stat'] = stat
                testing_res[method+'_pval'] = pval
            pvalue_res.append(testing_res)
            del testing_res
        ## data frame and calculate fdr
        pvalue_res = pd.DataFrame(pvalue_res)
        # pvalue_res = self._fdr_(pvalue_res = pvalue_res, testing_method = testing_method)
        return(pvalue_res)

    def _excu_commu_(self, s_m):
        """
        link function to self._commu_score_for_one_ and self._signif_test_ 
        """
        info('{}'.format(s_m))
        sensor, met = s_m.split(' ~ ')
        sensor_loc = self.avg_exp_indexer.tolist().index(sensor)
        avg_exp_sensor = self.avg_exp[sensor_loc].toarray()[0]
        
        met_loc = self.avg_met_indexer.tolist().index(met)
        avg_met_met = self.avg_met[met_loc].toarray()[0]
        
        sensor_perm_loc = self.permute_exp_indexer.tolist().index(sensor)
        perm_exp_avg_sensor = {i:self.perm_exp_avg[i][sensor_perm_loc].toarray()[0] for i in self.perm_exp_avg}
        
        met_perm_loc = self.permute_exp_indexer.tolist().index(sensor)
        perm_met_avg_met = {i:self.perm_met_avg[i][met_perm_loc].toarray()[0] for i in self.perm_met_avg}
        
        commu_mat, background_df = self._commu_score_for_one_(
                                    avg_exp_sensor = avg_exp_sensor, 
                                    avg_met_met = avg_met_met, 
                                    perm_exp_avg_sensor = perm_exp_avg_sensor, 
                                    perm_met_avg_met = perm_met_avg_met
                                   )
        pvalue_res = self._signif_test_(commu_mat = commu_mat, background_df = background_df)
        pvalue_res['Metabolite'] = met
        pvalue_res['Sensor'] = sensor
        
        return met+'~'+sensor, pvalue_res, background_df
    
                                   
    def _avg_exp_group_(self, group_names):
        ## avg exp by cell_group for met sensor
        avg_exp = np.empty(shape = (self.exp_mat.shape[0],0)) ## save exp data

        for x in group_names:
            cells = self.cell_ann[self.cell_ann['cell_group'] == x].index.tolist()
            cell_loc = [i for i, c in enumerate(self.exp_mat_columns) if c in cells]
            avg_exp = np.concatenate((avg_exp, self.exp_mat[:,cell_loc].mean(axis = 1)), axis = 1)
        self.avg_exp = sparse.csc_matrix(avg_exp)
        self.avg_exp_indexer = np.array(self.exp_mat_indexer)
        self.avg_exp_columns = np.array(group_names)
    
    def _avg_met_group_(self, group_names):
        """
        take average of sensor expression and metabolite by cell groups
        """
        ## avg met by cell_group for met
        avg_met = np.empty(shape = (self.met_mat.shape[0],0)) ## save exp data

        for x in group_names:
            cells = self.cell_ann[self.cell_ann['cell_group'] == x].index.tolist()
            cell_loc = [i for i, c in enumerate(self.met_mat_columns) if c in cells]
            avg_met = np.concatenate((avg_met, self.met_mat[:,cell_loc].mean(axis = 1)), axis = 1)
        self.avg_met = sparse.csc_matrix(avg_met)
        self.avg_met_indexer = np.array(self.met_mat_indexer)
        self.avg_met_columns = group_names

        
    def pred(self, n_shuffle = 1000, seed = 12345):
        """
        handling commmunication score
        """
        self.n_shuffle = n_shuffle
        self.seed = seed
        ## multiprocessing
        if self.thread is None:
            self.thread = 1
            
        info('Parameters: {shuffling: %s times, random seed: %s, thread: %s}'%(self.n_shuffle, self.seed, self.thread))
        
        ## avg exp and met level in cell group for real data
        group_names = self.cell_ann['cell_group'].unique().tolist()
        self.group_names = group_names
        
        ## make sure the same order by giving group_names
        if self.avg_met is None or self.avg_met_columns is not self.group_names:
            self._avg_met_group_(group_names = self.group_names)
        if self.avg_exp is None or self.avg_exp_columns is not self.group_names:
            self._avg_exp_group_(group_names = self.group_names)
            
        ## report data shape
        info('met_sensor: (%s, %s)'%self.met_sensor.shape)
        info('avg_exp: (%s, %s) for (gene, cell) of needed'%self.avg_exp.shape)
        info('avg_met: (%s, %s) for (metabolite, cell) of needed'%self.avg_met.shape)

        ## do n times shuffling for generating backgroud
        ## return a dict, index will be time of shuffling, given each shuffling, it is a matrix, rows are cell_group, columns are genes
        info('shuffling %s times for generating backgroud' % self.n_shuffle)
        self._shuffling_(group_names=self.group_names, n_shuffle = self.n_shuffle, seed = self.seed) 
        
        ## omputing commu score  for each met and sensor
        commu_res = collections.defaultdict() ## collect communication scores for plot later
        commu_res_bg = collections.defaultdict()

        info('thread: %s'%self.thread)
        pool = multiprocessing.Pool(self.thread)
        ## all sensor and met pairs
        s_m_all = self.met_sensor['Gene_name'] + ' ~ ' + self.met_sensor['HMDB_ID']
        res_col = pool.map(self._excu_commu_, s_m_all.tolist())
        pool.close()
        ## into dict
        for k, c, b in res_col:
            ## k is met~sensor, c is commu, b is background
            commu_res[k] = c
            commu_res_bg[k] = b
        ## dataframe of commu_res
        commu_res_df = pd.DataFrame()
        for pair in commu_res:
            commu_res_df = pd.concat([commu_res_df, commu_res[pair]])
        ## do fdr correction
        testing_method = ['ztest', 'ttest', 'ranksum_test', 'permutation_test']
        commu_res_df = self._fdr_(pvalue_res = commu_res_df, testing_method = testing_method)
        commu_res_df = commu_res_df.sort_values('permutation_test_fdr')
        return commu_res_df, commu_res_bg ## dict, keys are pair of met and sensor

        





