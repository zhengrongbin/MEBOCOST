#!/usr/bin/env python

# ================================
# @auther: Rongbin Zheng
# @email: Rongbin.Zheng@childrens.harvard.edu
# @date: May 2022
# ================================

import os,sys
import pandas as pd
import numpy as np
import collections
import pickle as pk
import traceback
from datetime import datetime
import xlmhg
from scipy import sparse

import multiprocessing

def info(string):
    """
    print information
    """
    today = datetime.today().strftime("%B %d, %Y")
    now = datetime.now().strftime("%H:%M:%S")
    current_time = today + ' ' + now
    print("[{}]: {}".format(current_time, string))


class PathwayEnrich:
    """
    class for pathway association to the communication
    briefly, the genes correlated with metabolite sensors and also high active in the cell might be the genes of associating with metabolite and communication

    Params
    -----
    commu_res
        a data frame of communication result.


    """
    def __init__(self,
                commu_res = pd.DataFrame,
                avg_exp = sparse.csc_matrix,
                avg_exp_indexer = np.array,
                avg_exp_columns = np.array,
                cell_ann = pd.DataFrame,
                gene_network = pd.DataFrame,
                gmt_path = None,
                min_term = 15,
                max_term = 500,
                thread = 2
                ):
        self.commu_res = commu_res
        self.avg_exp = avg_exp
        self.avg_exp_indexer = avg_exp_indexer
        self.avg_exp_columns = avg_exp_columns
        self.cell_ann = cell_ann
        self.gene_network = gene_network
        self.gmt_path = gmt_path
        self.min_term = min_term
        self.max_term = max_term
        self.thread = thread
        
#     def _avg_cell_(self):
#         """
#         take average of sensor expression and metabolite by cell groups
#         """
#         ## avg exp by cell_group for met sensor
#         info('Take average expression by cell group')
#         avg_exp = pd.merge(self.exp_mat, 
#                             self.cell_ann[['cell_group']], left_index = True, right_index = True).groupby('cell_group').mean()
#         return avg_exp


    def _read_gmt_(self):
        """
        read gmt file
        """
        info('Read gene set from GMT file')
        if self.gmt_path is None:
            raise KeyError('gmt_path is None')

        gmt = collections.defaultdict()
        with open(self.gmt_path) as f:
            for line in f:
                line = line.rstrip().split('\t')
                acc = line[0] ## pathway name or ID
                des = line[1] ## description of the pathway
                genes = line[2:] ## gene list
                if len(genes) <= self.max_term and len(genes) >= self.min_term:
                    gmt[acc] = genes 
        return(gmt)
        
    def _weight_gene_exp_(self):
        """
        weight the gene expression using gene network score 
        """
        info('Weight gene expression by gene network score of sensor')
        
        ## 1. gene cluster mean expression centerized to median by gene across cell type, to highlight the cell-type specificity
        avg_exp_norm = sparse.csr_matrix(list(map(lambda x: np.nan_to_num(x - np.median(x)), self.avg_exp.toarray())))
#         en = self.avg_exp.apply(lambda col: col - np.median(col)).T.dropna().T
        # 1.2 re-centerized to median by cell type across genes
        avg_exp_norm = sparse.csr_matrix(list(map(lambda x: np.nan_to_num(x - np.median(x)), avg_exp_norm.T.toarray()))).T
        self.avg_exp_norm = avg_exp_norm
        self.avg_exp_norm_indexer = self.avg_exp_indexer
        self.avg_exp_norm_columns = self.avg_exp_columns

#         en = en.apply(lambda row: row - np.median(row), axis = 1)
        ## 2. multiply weight 
        cg = list(set(self.gene_network.index.tolist()) & set(self.avg_exp_indexer.tolist())) ## common genes between gene ntk and exp
        weight = self.gene_network.reindex(index = cg)
        ## normed exp
        cg_loc = np.where(pd.Series(self.avg_exp_indexer).isin(cg))[0]
        cg_new = self.avg_exp_indexer[cg_loc]
        en = avg_exp_norm[cg_loc]
        en_indexer = cg_new
        en_columns = self.avg_exp_columns
        ## weight exp by gene network of each sensor
        weighted_exp = collections.defaultdict()
        ## function
        def _cal_func(e, ww, cutoff=0):
        #     e = (abs(e)/e) * np.sqrt(abs(e))
            e = np.nan_to_num(e)
            res = pd.Series()
            ## e>0, w>0
            index = (e>cutoff) & (ww>cutoff)
            res = pd.concat([res, (e[index] + ww[index])])
            ## e=0, w=0
            index = (e==0) | (ww==0)
            res = pd.concat([res, (e[index] * ww[index])])
            ## e<0, w<0
            index = (e<cutoff) & (ww<cutoff)
            res = pd.concat([res, (e[index] + ww[index])])
            ## e>0, w<0
            index = (e>cutoff) & (ww<cutoff)
            res = pd.concat([res, (-e[index] + ww[index])])
            ## e<0, w>0
            index = (e<cutoff) & (ww>cutoff)
            res = pd.concat([res, (e[index] - ww[index])])
            return(res)
        
        ## iterate for each sensor
        for s in self.commu_res['Sensor'].unique().tolist():
            w = weight[s]
            w = w[~np.isinf(w)]
            genes = w.index.tolist()
            ## finial weighted exp
            en_gene_loc = np.where(pd.Series(en_indexer).isin(genes))[0]
            en_genes = en_indexer[en_gene_loc]
            exp = en[en_gene_loc].T ## rows = cell, columns = gene
            ## reorder w
            w = w[en_genes]
            ## for each cell type, get weight for all genes
            #.apply(lambda row: _cal_func(row, w), axis = 1)
            w_exp = [_cal_func(l.toarray()[0], w) for l in exp]
            w_exp = pd.DataFrame(w_exp, index = en_columns)
             
            weighted_exp[s] = {'weight_exp': sparse.csr_matrix(w_exp),
                              'index': np.array(w_exp.index),
                              'columns': np.array(w_exp.columns)}
            
        return(weighted_exp)


    def getmHG(self, genes_indices, N, X, L):
        # genes_indices = sorted(set(genes_indices))
        go_indices = np.uint16(genes_indices)
        res = xlmhg.get_xlmhg_test_result(N, go_indices, X=X, L=L)
        try:
            res_list = ["%s,%s,%s,%s"%(res.N, res.cutoff, res.K, res.k), res.stat, res.fold_enrichment, res.pval]
        except:
            res_list = ["%s,%s,%s,%s"%(res.N, res.cutoff, res.K, res.k), res.stat, 0, res.pval]
        # res_list = [res.N, res.cutoff, res.K, res.k, res.fold_enrichment, res.pval]
        return [res]+ res_list
    
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
    
    ## the maximum length of a term is 374 genes, therefore, roughly set L to 500 is enough
    def _link_mHG(self, weight, gmt, dL):
        w = weight.sort_values(ascending = False).copy()
        genes = []
        for x in gmt:
            genes.extend(gmt[x])
        genes = list(set(genes))

        indices = {x:np.where(w.index.isin(gmt[x]))[0] for x in gmt}
        N = w.size
        L = min(dL, len(set(genes) & set(w.index.tolist()))) ## dL is the default int for upper limit of L
        res = pd.DataFrame(map(lambda x: self.getmHG(indices[x], N=int(N), X = 1, L = L), indices),
                      index = list(indices.keys()))
        res.columns = ['obj', 'N,B,n,b', 'Stats', 'FoldEnrichment', 'pval']
        res['fdr'] = self.bh_fdr(res['pval'])
        resobjs = dict(res['obj'])
        res = res.drop('obj', axis = 1)
        res = res.sort_values('fdr')
        ## length
        res['gsLength'] = [len(set(gmt[x]) & set(w.index.tolist())) for x in res.index.tolist()]
        return(resobjs, res)
    

    def _excu_sensor_enrich_(self, s_r):
        """
        link enrichment for each
        """
        info(s_r)
        # print(self.weighted_exp)
        s, r = s_r.split(' ~ ') ## sensor ~ receiver
        
        s_loc = np.where(self.weighted_exp[s]['index'] == r)
        w_e = pd.Series(self.weighted_exp[s]['weight_exp'][s_loc].toarray()[0],
                               index = self.weighted_exp[s]['columns']).sort_values(ascending = False) ## weight exp for each pair of sensor and receiver
        if s in w_e.index.tolist():
            w_e = w_e.drop(s)
        w_e = w_e[~np.isinf(w_e)].sort_values(ascending = False)
#         print(w_e)
#         print(self.gmt)
        resobjs, res = self._link_mHG(w_e, self.gmt, 500)
#         resobjs, res = self._link_mHG(weight = w_e, gmt = self.gmt, dL = 500)
        res = res.sort_values('fdr')
        return(s_r, resobjs, res)

    def _filter_focus_(self):
        """
        focus on some communications if given
        """
        focus_commu = self.commu_res.copy()
        if self.sender_focus:
            focus_commu = focus_commu[focus_commu['Sendor'].isin(self.sender_focus)]
        if self.metabolite_focus:
            focus_commu = focus_commu[focus_commu['Metabolite_Name'].isin(self.metabolite_focus) | focus_commu['Metabolite'].isin(self.metabolite_focus)]
        if self.sensor_focus:
            focus_commu = focus_commu[focus_commu['Sensor'].isin(self.sensor_focus)]
        if self.receiver_focus:
            focus_commu = focus_commu[focus_commu['Reciver'].isin(self.receiver_focus)]
        return focus_commu
        
    def _sensor_enrich_(self):
        """
        weight_exp
            a pd.Series, index = gene name, value = weighted expression
        """
        info('Enrichment for significant sensor in receiver cell')
        comm_df = self._filter_focus_()
        if comm_df.shape[0] == 0:
            info('No Sensor in the focus')
        sr = np.unique(comm_df['Sensor'] + ' ~ ' + comm_df['Receiver']) ## all sensor in receivers
        info('Thread: %s'%(self.thread))
        pool = multiprocessing.Pool(self.thread)
        res_col = pool.map(self._excu_sensor_enrich_, sr)
        pool.close()
        ## collect
        sensor_enrich_res = {s_r: {'mHG_obj':resobjs, 'mHG_res':res} for s_r, resobjs, res in res_col}

        return(sensor_enrich_res)

    
    def _deconv_(self, s_r):
        """
        deconvolution
        """
        s, r = s_r.split(' ~ ')
        ## weight exp for each pair of sensor and receiver
        s_loc = np.where(self.weighted_exp[s]['index'] == r)
        w_e = pd.Series(self.weighted_exp[s]['weight_exp'][s_loc].toarray()[0],
                               index = self.weighted_exp[s]['columns'])
        
        if s in w_e.index.tolist():
            w_e = w_e.drop(s)
        w_e = w_e[~np.isinf(w_e)]
        ## deconvolute communications
        comm_w = self.commu_res.loc[(self.commu_res['Sensor'] == s) &
                                    (self.commu_res['Receiver'] == r),'commu_weight']

        comm_frac = self.comm_frac_w.loc[(self.comm_frac_w['Sensor'] == s) &
                                  (self.comm_frac_w['Receiver'] == r), 'frac']
        ## pairwise multiplication gene weight * cell pair
        deconv = pd.DataFrame([[i * j for j in comm_frac] for i in w_e],
                              index = w_e.index, columns = comm_frac.index).T
        deconv.index = deconv.index + ' ~ ' + s_r ## x is sensor and receiver pair
        return deconv

    def _deconvolution_and_agg_(self):
        """
        deconvolution of sensor weighted expression by communication events
        the assumption is that the activation of sensor-related genes are regulated by all communications to the corresponding sensor
        """
        info('Weighted expression deconvolution to metabolite-sensor events')
        ## deconvolution of weighted exp into communications (sensor-receiver same, but senders different, actually  (sender - metabolite - sensor - receiver))
        sensor_receptors = np.unique(self.commu_res['Sensor'] + ' ~ ' + self.commu_res['Receiver'])
        info('{} sensor-receiver pairs'.format(len(sensor_receptors)))
        
        ## commu weight as score * -log10(p)
        self.commu_res['commu_weight'] = self.commu_res['Commu_Score'] #-np.log10(self.commu_res[self.pval_method]) * self.commu_res['Commu_Score']
        self.commu_res.index = range(self.commu_res.shape[0])
        
        ## the weight of save sensor divided to each commu (sender-met)
        comm_frac_w = self.commu_res.groupby(['Sensor', 'Receiver']).apply(lambda df: df['commu_weight'] / df['commu_weight'].sum()).reset_index()
        comm_frac_w.columns = ['Sensor','Receiver','index','frac']
        ## cat more info
        comm_frac_w['Metabolite_Name'] = self.commu_res.reindex(comm_frac_w['index'])['Metabolite_Name'].tolist()
        comm_frac_w['Sender'] = self.commu_res.reindex(comm_frac_w['index'])['Sender'].tolist()
        comm_frac_w.index = comm_frac_w['Sender'] + ' ~ ' + comm_frac_w['Metabolite_Name']
        self.comm_frac_w = comm_frac_w

        pool = multiprocessing.Pool(self.thread)
        res_col = pool.map(self._deconv_, sensor_receptors)
        pool.close()
        ## concat to data frame
        w_e_d = pd.DataFrame()
        for d in res_col:
            w_e_d = pd.concat([w_e_d, d])

        info('Weighted expression deconvolution to sender-receiver events')
        ### aggerate weighted exp by sender-receiver pair
        w_e_a = collections.defaultdict()
        sr_pair = np.unique(self.commu_res['Sender'] + ' ~ ' + self.commu_res['Receiver'])
        info('{} sender-receiver pairs'.format(len(sr_pair)))
        for pair in sr_pair:
            sender, receiver = pair.split(' ~ ')
            ## mean weighted exp in multiple communication
            tmp_exp_index = w_e_d.index.str.startswith(sender) & w_e_d.index.str.endswith(receiver)
            ## take the mean for multiple comunications for a same cell pair (sender and receiver)
            tmp_exp = w_e_d.loc[tmp_exp_index].mean().sort_values(ascending = False)
            w_e_a[pair] = tmp_exp
        w_e_a = pd.DataFrame(w_e_a)

        return w_e_d, w_e_a


    def _excu_cell_enrich_(self, p):
        ## p is the cell pair
        info(p)
        geneScore = self.weighted_exp_agg[p].dropna().sort_values(ascending = False)
        resobjs, res = self._link_mHG(weight = geneScore, gmt = self.gmt, dL = 500)
        res = res.sort_values('fdr')
        return(p, resobjs, res)


    def _cell_enrich_(self):
        """
        pathway enrichment for a pair of cell-cell communication
        """
        info('Enrichment for cell-cell communication events')
        ## all cell pairs
        comm_df = self._filter_focus_()
        if comm_df.shape[0] == 0:
            info('No cells in the focus')
        
        cr = np.intersect1d(self.weighted_exp_agg.columns, np.unique(comm_df['Sender']+' ~ '+comm_df['Receiver']))
        info('Thread: %s'%(self.thread))
        pool = multiprocessing.Pool(self.thread)
        res_col = pool.map(self._excu_cell_enrich_, list(cr))
        pool.close()
        ## collect
        cell_enrich_res = {c_r: {'mHG_obj':resobjs, 'mHG_res':res} for c_r, resobjs, res in res_col}

        return cell_enrich_res


    def _load_data_(self, pval_method = 'ranksum_test'):
        """
        load data and calculate weighted expression
        """
        self.pval_method = pval_method
#         ## take average gene expression across cell group
#         self.avg_exp = self._avg_cell_() ## cell_group by gene matrix
        ## read gmt
        self.gmt = self._read_gmt_()
        ## select terms
        self.gmt = {x: self.gmt[x] for x in self.gmt.keys() if len(self.gmt[x]) > self.min_term and len(self.gmt[x]) < self.max_term for x in self.gmt.keys()}
        ## weight gene exp by gene network 
        self.weighted_exp = self._weight_gene_exp_()
        ## deconvolution of sensor-gene associtation into communication events (sender - metabolite - sensor - receiver)
        self.weighted_exp_deconv, self.weighted_exp_agg = self._deconvolution_and_agg_()
        return

    def _clean_return_dict_(self):
        """
        return the dict as result
        """
        res = collections.defaultdict()

        res['parameters'] = {'gmt_path':self.gmt_path,
                            'min_term':self.min_term,
                            'max_term':self.max_term,
                            'thread':self.thread,
                            'pval_method':self.pval_method,
                            'sender_focus':self.sender_focus,
                            'metabolite_focus':self.metabolite_focus,
                            'sensor_focus':self.sensor_focus,
                            'receiver_focus':self.receiver_focus}
        res['gmt'] = self.gmt
        res['avg_exp_norm'] = pd.DataFrame(self.avg_exp_norm.toarray(),
                                          index = self.avg_exp_norm_indexer,
                                          columns = self.avg_exp_norm_columns)
#         res['avg_exp_norm_indexer'] = self.avg_exp_norm_indexer
#         res['avg_exp_norm_columns'] = self.avg_exp_norm_columns

        res['weighted_exp'] = self.weighted_exp
        res['comm_frac_w'] = self.comm_frac_w
        res['weighted_exp_deconv'] = self.weighted_exp_deconv
        res['weighted_exp_agg'] = self.weighted_exp_agg
        res['sensor_res'] = self.sensor_res
        res['cellpair_res'] = self.cellpair_res
        return res
        
        
    
    def _pred_(self, pval_method = 'ranksum_test', 
               sensor_in_receiver = True, 
               sender_to_receiver = True, 
               Return = False,
               sender_focus = [],
               metabolite_focus = [],
               sensor_focus = [],
               receiver_focus = []):
        """
        infer pathways of those are associated with communications either sensor or cell-cell
        """
        ### load data
        self._load_data_(pval_method = pval_method)
        ## filter
        self.sender_focus = sender_focus
        self.metabolite_focus = metabolite_focus
        self.sensor_focus = sensor_focus
        self.receiver_focus = receiver_focus
        ## pathway enrichment for significant sensor in receiver
        sensor_res = collections.defaultdict()
        if sensor_in_receiver:
            ## iterate for each sensor in receivers
            sensor_res = self._sensor_enrich_()
        
        ## pathway enrichment for significant sender-receiver pairs
        cellpair_res = collections.defaultdict()
        if sender_to_receiver:
            ## pathway enrichment for cell-cell pair
            cellpair_res = self._cell_enrich_()

        self.sensor_res = sensor_res
        self.cellpair_res = cellpair_res
        
        if Return:
            return self._clean_return_dict_()






