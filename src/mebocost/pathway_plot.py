#!/usr/bin/env python

# ================================
# @auther: Rongbin Zheng
# @email: Rongbin.Zheng@childrens.harvard.edu
# @date: May 2022
# ================================

import os,sys
import time
from datetime import datetime
import numpy as np
import pandas as pd
import traceback
import matplotlib
import collections
import xlmhg
import seaborn as sns
from adjustText import adjust_text
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib.patches import Rectangle
import networkx as nx
import re
## disable warnings
import warnings
warnings.filterwarnings("ignore")


plt.rcParams.update(plt.rcParamsDefault)
rc={"axes.labelsize": 16, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.titleweight":"bold", #"font.size":14,
    "figure.figsize":(5.5,4.2), "font.weight":"regular", "legend.fontsize":10,
    'axes.labelpad':8, 'figure.dpi':300}
plt.rcParams.update(**rc)


def info(string):
    """
    print information
    """
    today = datetime.today().strftime("%B %d, %Y")
    now = datetime.now().strftime("%H:%M:%S")
    current_time = today + ' ' + now
    print("[{}]: {}".format(current_time, string))
    return

def _select_pathway_(res_dict,
                    a_pair,
                    pval_cutoff = 0.05, 
                    ES_cutoff = 0,
                    maxSize = 500,
                    minSize = 15):
    """
    select pathway based on cutoff, and make to dataframe
    """
    plot_df = res_dict[a_pair]['mHG_res']
    if ES_cutoff >= 0:
        plot_df = plot_df[(plot_df['fdr'] < pval_cutoff) &
                         (plot_df['FoldEnrichment'] > ES_cutoff) &
                         (plot_df['gsLength'] >= minSize) &
                         (plot_df['gsLength'] <= maxSize)]
    else:
        plot_df = plot_df[(plot_df['fdr'] < pval_cutoff) &
                         (plot_df['FoldEnrichment'] < ES_cutoff) &
                         (plot_df['gsLength'] >= minSize) &
                         (plot_df['gsLength'] <= maxSize)]
    
    plot_df['pathway'] = [x.replace(re.search('hsa[0-9]* |mmu[0-9]* ', x).group(), '') for x in plot_df.index.tolist()]
#     ## remove ambiguous, namely predicted as both positive ane negative 
#     duplicated_pathways = plot_df['pathway'][plot_df['pathway'].duplicated()]
#     plot_df = plot_df[~plot_df['pathway'].isin(duplicated_pathways)]
    ## sorted by ES
    plot_df = plot_df.sort_values('FoldEnrichment')
    return(plot_df)



def _scatter_(res_dict,
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
                pdf = None,
                show_plot = False,
                return_fig = False):
    """
    This function to draw a scatter plot to show significant pathways that strongly associated with a certain communication, can be sensor in receiver or sender to receiver
    
    Params
    ------
    result_dict: a dict, containing communication result, usually it is either gsea_res_cell_cell or gsea_res in config output
    a_pair: the key of result_dict, should be such as Malignant_4->Malignant_0 in gsea_res_cell_cell and SLC38A1~CD8T_6 in gsea_res
    """
    ##clean
    plt.close() 
    ## plot for one enrichment
    # cell_pair = 'Malignant_4->Malignant_0'
    plot_df = _select_pathway_(
                            res_dict=res_dict,
                            a_pair=a_pair,
                            pval_cutoff = pval_cutoff, 
                            ES_cutoff = ES_cutoff,
                            maxSize = maxSize,
                            minSize = minSize
                            )
    if plot_df.shape[0] == 0:
        info('No engouth pathway to plot!')
        return
        
    plot_df['fdr'] = plot_df['fdr'].astype('float')
    ## if the scale if too small
    logp = -np.log10(plot_df['fdr'])
    if (max(logp) - min(logp)) < 0.1:
        vmin = -np.log10(0.05)
    
    if plot_df.shape[0] == 0:
        info('no enough significant pathway to show!')
        return
    if not vmax:
        vmax = np.percentile(-np.log10(plot_df['fdr']), 75)
    if not vmin:
        vmin = np.percentile(-np.log10(plot_df['fdr']), 25)

    nrow, ncol = plot_df.shape
    if figsize == 'auto':
        figsize = (5+nrow*0.25, 3+nrow*0.16)

    fig, ax = plt.subplots(figsize = figsize)
    sp = ax.scatter(x = plot_df['FoldEnrichment'],
                    y = plot_df['pathway'],
                    s = plot_df['gsLength'],
                    c = -np.log10(plot_df['fdr']),
                    cmap = cmap, zorder = 100,
                   vmax = vmax, vmin = vmin,
                   edgecolor = 'none',
                   alpha = .8)
    cbar = plt.colorbar(sp, ax = ax, 
                        shrink = .4,
                        location = 'right')
    cbar.set_label('-log10(p.adj)', fontsize = 10)
    ## size legend
    slegend = [np.percentile(plot_df['gsLength'], 25),
              np.percentile(plot_df['gsLength'], 50),
              np.percentile(plot_df['gsLength'], 75)]
    for i in slegend:
        ax.scatter([], [],
                  s = int(i),
                   facecolor = 'none',
                   edgecolor = 'black',
                  label = int(i))
    ax.legend(title = 'Gene Count',
                  loc = 'upper left',
                 bbox_to_anchor=(1.02, 1),
                 frameon = False, 
              fontsize = 10, title_fontsize = 10)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_ylim(-0.5, nrow-.5)
    ax.set_xlabel('Fold Enrichment')
    ax.set_title(title, fontsize = 14, pad = 20, fontweight = 'bold')
    ax.grid()
    plt.tight_layout()
    pdf.savefig(fig) if pdf else None
    if show_plot:
        plt.show()
    plt.close()
    if return_fig:
        return(fig)
    

def _stacked_bar_(res_dict,
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
                pdf = None,
                show_plot = False,
                return_fig = False):
    """
    this function to draw figure to compare potential pathways between two conditions, the conditions can be two pairs of cells
    -------------
    result_dict: a dict, containing communication result, usually it is either gsea_res_cell_cell or gsea_res in config output
    pair1: the key of result_dict, should be such as Malignant_4->Malignant_0 in gsea_res_cell_cell and SLC38A1~CD8T_6 in gsea_res
    pair2: same format with pair1
    """
    ##clean
    plt.close()
    plot_df1 = _select_pathway_(
                            res_dict=res_dict,
                            a_pair=pair1,
                            pval_cutoff = pval_cutoff, 
                            ES_cutoff = ES_cutoff,
                            maxSize = maxSize,
                            minSize = minSize
                            )

    if plot_df1.shape[0] == 0:
        info('no significant enrichment for pair, please try one pair plot by calling _scatter_ function')
        return
            
    plot_df2 = _select_pathway_(
                            res_dict=res_dict,
                            a_pair=pair2,
                            pval_cutoff = pval_cutoff, 
                            ES_cutoff = ES_cutoff,
                            maxSize = maxSize,
                            minSize = minSize
                            )
    if plot_df2.shape[0] == 0:
        info('no significant enrichment for pair, please try one pair plot by calling _scatter_ function')
        return
    
    ## merge two df for plot
    plot_df = pd.merge(plot_df1[['pathway', 'FoldEnrichment', 'fdr']], 
                         plot_df2[['pathway', 'FoldEnrichment', 'fdr']], 
                         left_on = 'pathway',
                         right_on = 'pathway', how = 'outer').fillna(0)
    plot_df.columns = ['pathway', 'cp1_FE', 'cp1_fdr', 'cp2_FE', 'cp2_fdr']
#     print(plot_df)
    if plot_df.shape[0] == 0:
        info('no enough significant pathway to compare')
        return

    plot_df = pd.concat([
        plot_df[plot_df['cp1_FE'] == 0].sort_values('cp2_FE'),
        plot_df[(plot_df['cp1_FE'] != 0) & (plot_df['cp2_FE'] != 0)].sort_values('cp2_FE'),
        plot_df[plot_df['cp2_FE'] == 0].sort_values('cp1_FE')
    ])
    #### ==== bar plot ======
    if not vmax:
        vmax = np.percentile([-np.log10(p) for p in plot_df['cp1_fdr'].tolist()+plot_df['cp2_fdr'].tolist() if p != 0], 75)
    if not vmin:
        vmin = np.percentile([-np.log10(p) for p in plot_df['cp1_fdr'].tolist()+plot_df['cp2_fdr'].tolist() if p != 0], 25)

    nrow, ncol = plot_df.shape
    if figsize == 'auto':
        figsize = (8+nrow*0.2, 3.5+nrow*0.17)

    fig, ax = plt.subplots(constrained_layout=True, ncols=2, nrows = 1, 
                           sharey = True, 
                           gridspec_kw = {'width_ratios':[1,1],
                                         'wspace':.01},
                           figsize = figsize)
    ## left side
    ssum = plot_df[['cp1_FE', 'cp2_FE']].T.sum()
    stacked_df = plot_df[['cp1_FE', 'cp2_FE']].apply(lambda col: col / ssum)
    stacked_df.index = plot_df['pathway'].tolist()
    yticks = list(np.arange(0, nrow * 3, 3) + 0.5)

    ax[0].barh(yticks, stacked_df['cp1_FE'], 
               height = 1.2, label = pair1,
              color = colors[0])
    ax[0].barh(yticks, stacked_df['cp2_FE'],
               height = 1.2,
               left=stacked_df['cp1_FE'],
               label = pair2,
               color = colors[1])
    ax[0].vlines(0.5, *ax[0].get_ylim(), linestyle = 'dashed', color = 'lightgrey')
    ax[0].legend()
    ax[0].set_xlabel('Relative Fold Enrichment')

    # plt.yticks(ticks = yticks, labels=stacked_df.index)


    ## right side
    my_cmap = plt.cm.get_cmap(cmap)
    norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
    barcolor1 = np.array(list(map(norm, [-np.log10(p) if p != 0 else 0 for p in plot_df['cp1_fdr'].tolist()])))
    barcolor1 = [my_cmap(x) for x in barcolor1]
    ax[1].barh(np.arange(0, nrow * 3, 3)-0.05, 
            plot_df['cp1_FE'], 
            color = barcolor1,
           edgecolor = 'lightgrey',
           label = pair1)

    barcolor2 = np.array(list(map(norm, [-np.log10(p) if p != 0 else 0 for p in plot_df['cp2_fdr'].tolist()])))
    barcolor2 = [my_cmap(x) for x in barcolor2]
    ax[1].barh(np.arange(1, nrow * 3, 3)+0.05, 
            plot_df['cp2_FE'],
           color = barcolor2,
           edgecolor = 'black',
           label = pair2)

    plt.yticks(ticks = yticks,
                 labels = plot_df['pathway'].tolist(), size = 12)

    sm =matplotlib.cm.ScalarMappable(cmap=my_cmap,norm=norm)
    sm.set_array([])
    cbar=plt.colorbar(sm,label='',shrink=.5)
    cbar.set_label('-log10(p.adj)', fontsize = 10)
    ## bar group legend
    ax[1].legend(bbox_to_anchor=(1, .9),
             frameon = False, fontsize = 10)
    ax[1].set_ylim(-1, nrow*3+1)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].set_xlabel('Fold Enrichment')
    ax[1].set_title(title, fontsize = 12, pad = 20, fontweight = 'bold')  if title != '' else None

    # plt.tight_layout()
    pdf.savefig(fig) if pdf else None
    if show_plot:
        plt.show()
    plt.close()
    if return_fig:
        return(fig)


def _multi_dot_(res_dict,
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
                    pdf = None,
                    show_plot = False,
                    swap_axis = False,
                    return_fig = False):
    """
    this function to draw figure to compare potential pathways in multiple conditions, the conditions can be many pairs of cells
    -------------
    result_dict: a dict, containing communication result, usually it is either gsea_res_cell_cell or gsea_res in config output
    pairs: a list, cell pairs in format such as ['CD8T_6->Malignant_0', 'Malignant_4->Malignant_0', 'Endothelial_10->Malignant_0', 'CD4Tconv_15->Malignant_0']
    flip: flip x and y axis
    """    
    ##clean
    plt.close() 
    plot_df = pd.DataFrame()
    for cp in pairs:
        tmp_df = _select_pathway_(
                            res_dict=res_dict,
                            a_pair=cp,
                            pval_cutoff = pval_cutoff, 
                            ES_cutoff = ES_cutoff,
                            maxSize = maxSize,
                            minSize = minSize
                            )
        tmp_df['pair'] = cp
        if tmp_df.shape[0] == 0:
            info('no significant enrichment for {}'.format(cp))
        plot_df = pd.concat([plot_df, tmp_df])
    if plot_df.shape[0] == 0:
        info('no enough significant pathway to plot')
        return

    ## order pathway
    pathway_order = []
    for pathway in plot_df['pathway'].unique():
        cps = plot_df[plot_df['pathway']==pathway]['pair'].unique().tolist()
        pathway_order.append([pathway, ';'.join(sorted(cps)), len(cps)])
    pathway_order = pd.DataFrame(pathway_order)
    ## order for cell pair
    cp_order = pathway_order[pathway_order[2]==1].groupby(1)[1].count().reindex(pairs).fillna(0).sort_values(ascending = False).index.tolist()
    pathway_order = pathway_order.sort_values([2,1])[0].tolist()

    ## force the order
    plot_df['pathway'] = plot_df['pathway'].astype('category').cat.set_categories(pathway_order)
    plot_df.sort_values(by=['pathway'], inplace=True)
    plot_df['pair'] = plot_df['pair'].astype('category').cat.set_categories(cp_order)
    plot_df.sort_values(by=['pathway'], inplace=True)

    nrow = len(pathway_order)
    ncol = len(plot_df['pair'].unique())

    if figsize == 'auto':
        figsize = (5+ncol*0.25, 7+nrow*0.1) if not swap_axis else (7+nrow*0.15, 5+ncol*0.25)

    node_size_norm_fun = lambda x, y: node_size_norm[0]+((x-min(y)) / (max(y) - min(y)) * (node_size_norm[1]-node_size_norm[0])) if max(y) != min(y) else node_size_norm[0]+((x-min(y)) / max(y) * (node_size_norm[1]-node_size_norm[0])) 

    if not vmax:
        vmax = np.percentile(-np.log10(plot_df['fdr']), 75)
    if not vmin:
        vmin = np.percentile(-np.log10(plot_df['fdr']), 25)

    fig, ax = plt.subplots(figsize = figsize)

    sp = ax.scatter(x = plot_df['pair'] if not swap_axis else plot_df['pathway'],
              y = plot_df['pathway'] if not swap_axis else plot_df['pair'],
              marker = 's',
              c = -np.log10(plot_df['fdr']),
              s = list(map(lambda x: node_size_norm_fun(x, plot_df['FoldEnrichment']), plot_df['FoldEnrichment'])) if plot_df['FoldEnrichment'].tolist()[0] > 0 else list(map(lambda x: node_size_norm_fun(abs(x), abs(plot_df['FoldEnrichment'])), abs(plot_df['FoldEnrichment']))),
              cmap = cmap,
              vmax = vmax,
              vmin = vmin,
              alpha = .8,
              zorder = 100)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_ylim(-0.5, nrow-.5) if not swap_axis else ax.set_ylim(-0.5, ncol-.5)
    ax.set_xlim(-0.5, ncol-.5) if not swap_axis else ax.set_xlim(-0.5, nrow-.5)
    ax.grid(color = 'lightgrey')
    ax.tick_params(axis = 'x', rotation = 90)
    ax.tick_params(axis = 'y', labelsize = 10)

    cbar = plt.colorbar(sp, ax = ax, 
                        shrink = .4 if not swap_axis else .8,
                        location = 'right')
    cbar.set_label('-log10(p.adj)', fontsize = 10)

    ## size legend
    slegend = [np.percentile(plot_df['FoldEnrichment'], 1),
              np.percentile(plot_df['FoldEnrichment'], 50),
              np.percentile(plot_df['FoldEnrichment'], 100)]
    for i in slegend:
        ax.scatter([], [],
                   s = node_size_norm_fun(i, plot_df['FoldEnrichment']) if plot_df['FoldEnrichment'].tolist()[0] > 0 else node_size_norm_fun(abs(i), abs(plot_df['FoldEnrichment'])),
                   facecolor = 'none',
                   edgecolor = 'black',
                   marker = 's',
                   label = round(i, 2))
    ax.legend(title = 'Fold Enrichment',
              loc = 'upper left',
              bbox_to_anchor=(1.02, 1) if not swap_axis else (1.15, 1),
              frameon = False, 
              fontsize = 10,
              title_fontsize = 10)
    ax.set_title(title, fontsize = 12, pad = 20, fontweight = 'bold') if title != '' else None
    plt.tight_layout()
    pdf.savefig(fig) if pdf else None
    if show_plot:
        plt.show()
    plt.close()
    if return_fig:
        return(fig)

def is_equal(a, b, tol):
    """Ratio test to check if two floating point numbers are equal.
    Parameters
    ----------
    a: float
        The first floating point number.
    b: float
        The second floating point number.
    tol: float
        The tolerance used.
    Returns
    -------
    bool
        Whether or not the two numbers are deemed equal.
    """
    if a == b or abs(a-b) <= tol * max(abs(a), abs(b)):
        return True
    else:
        return False
    
def get_hgp(p, k, N, K, n):
    """Calculate the hypergeometric p-value when p = f(k; N,K,n) is already known.
    """
    pval = p
    while k < min(K, n):
        p *= (float((n-k)*(K-k) / float((k+1)*(N-K-n+k+1))))
        pval += p
        k += 1
    return pval

def get_hypergeometric_stats(N, indices):
    """Calculates hypergeom. p-values and fold enrichments for all cutoffs.
    Parameters
    ----------
    N: int
        The length of the list
    indices:  `numpy.ndarray` with ``dtype=np.uint16``
        The (sorted) indices of the "1's" in the list.
    """
    assert isinstance(N, (int, np.integer))
    assert isinstance(indices, np.ndarray) and \
            np.issubdtype(indices.dtype, np.uint16)

    K = indices.size

    pvals = np.empty(N+1, dtype=np.float64)
    folds = np.empty(N+1, dtype=np.float64)
    pvals[0] = 1.0
    folds[0] = 1.0

    n = 0
    k = 0
    p = 1.0
    while n < N:
        if k < K and indices[k] == n:
            # "add one"
            # calculate f(k+1; N,K,n+1) from f(k; N,K,n)
            p *= (float((n+1) * (K-k)) / \
                  float((N-n) * (k+1)))
            k += 1
        else:
            # "add zero"
            # calculate f(k; N,K,n+1) from f(k; N,K,n)
            p *= (float((n+1) * (N-K-n+k)) /
                  float((N-n) * (n-k+1)))
        n += 1
        # calculate hypergeometric p-value
        pvals[n] = get_hgp(p, k, N, K, n)
        # calculate fold enrichment
        folds[n] = k / (K*(n/float(N)))

    return pvals, folds

def _text_quadrant(expList, gxtList, geneSet):
    mx = np.median(expList)
    q1gene = expList[geneSet][(gxtList[geneSet]>0) & (expList[geneSet]>mx)]
    q2gene = expList[geneSet][(gxtList[geneSet]>0) & (expList[geneSet]<mx)]
    q3gene = expList[geneSet][(gxtList[geneSet]<0) & (expList[geneSet]<mx)]
    q4gene = expList[geneSet][(gxtList[geneSet]<0) & (expList[geneSet]>mx)]
    res=(len(q1gene), len(q2gene), len(q3gene), len(q4gene))
    return(res)
    
def _sensor_ES_plot_(geneSet, 
                  mHG_obj,
                  expList,
                  gxtList,
                  sensor = '',
                  receiver = '',
                  figsize = (8, 3.5),
                  dot_color = '#1874CD',
                  curve_color = 'black',
                  title = '',
                  pdf = None,
                  show_plot = False,
                  return_fig = False
             ):
    ##clean
    plt.close() 
    ## geneset
    cgene = np.intersect1d(expList.index, gxtList.index)
    geneSet = np.intersect1d(geneSet, cgene)
    ## mHG obj
    N = mHG_obj.N
    indices = mHG_obj.indices
    pvals, folds = get_hypergeometric_stats(mHG_obj.N, mHG_obj.indices)
    try:
        fold, pval = mHG_obj.fold_enrichment, mHG_obj.pval
    except:
        fold, pval = 0, mHG_obj.pval
    ## figure
    fig = plt.figure(constrained_layout=True, figsize = figsize)
    subfigs = fig.add_gridspec(2, 2, width_ratios=[1,1],
                               height_ratios = [1,0.1],
                               hspace = 0, wspace = 0.05)
    left = fig.add_subplot(subfigs[:, 0])
    right_upper = fig.add_subplot(subfigs[0, 1])
    right_lower = fig.add_subplot(subfigs[1, 1], sharex = right_upper)
    ## scatter
    ## left scatter
    mx = np.median(expList)
    left.scatter(x = expList[cgene], y = gxtList[cgene], s = 1, 
               edgecolor = 'none', facecolor = 'lightgrey')
    ymin, ymax = left.axes.get_ylim()
    left.vlines(mx, ymin, ymax, color = 'black',
              linestyle = '--', linewidth = .7)
    xmin, xmax = left.axes.get_xlim()
    left.hlines(0, xmin, xmax, color = 'black',
              linestyle = '--', linewidth = .7)
    left.scatter(x = expList[geneSet], y = gxtList[geneSet], 
               color = dot_color, s = 10, alpha = .8)
    ## ann
    q1, q2, q3, q4 = _text_quadrant(expList, gxtList, geneSet)
    
    extend = 0
    left.text(x = xmax, y = ymax - extend,
            s = '{}/{}'.format(q1,len(geneSet)),
            color = dot_color, ha = 'right')
    left.text(x = xmin, y = ymax - extend,
            s = '{}/{}'.format(q2,len(geneSet)),
            color = dot_color, ha = 'left')
    left.text(x = xmin, y = ymin + extend,
            s = '{}/{}'.format(q3,len(geneSet)),
            color = dot_color, ha = 'left')
    left.text(x = xmax, y = ymin + extend,
            s = '{}/{}'.format(q4,len(geneSet)),
            color = dot_color, ha = 'right')
    if receiver != '':
        left.set_xlabel('Scaled Expression in %s'%receiver, fontsize = 12, labelpad = 10)
    else:
        left.set_xlabel('Scaled Expression in Receiver', fontsize = 12, labelpad = 10)
    if sensor != '':
        left.set_ylabel('Gene Correlation with %s'%sensor, fontsize = 12)
    else:
        left.set_ylabel('Gene Correlation with Sensor', fontsize = 12)
#     sns.despine(trim = True, offset = 2, ax = left)
    left.spines['right'].set_visible(False)
    left.spines['top'].set_visible(False)
    ## enrich curve
    right_upper.plot(range(N+1), -np.log10(pvals), #marker = 1,
               color = curve_color, markersize = 5, linewidth = .7)
    right_upper.scatter(indices, -np.log10(pvals)[indices],
                        color = dot_color, s = 5, zorder = 3,
                      label = 'Pathway gene')

    right_upper.set_ylabel('Enrichment score', fontsize = 12)
    right_upper.spines['right'].set_visible(False)
    right_upper.spines['top'].set_visible(False)
#     right_upper.spines['bottom'].set_visible(False)
    right_upper.set_xticks([])
    ymin, ymax = right_upper.get_ylim()
    right_upper.text(0, ymin-(ymax*0.065), 'Highly associated', color = 'red')
    right_upper.text(N, ymin-(ymax*0.065), 'Lowly associated', 
               color = 'blue', ha = 'right')
    leg = right_upper.legend(fontsize = 10, frameon = False)
    right_upper.set_xlim(-100, N)
    ## label scores
    x = (N / 3) * 1.7
    y = (right_upper.get_ylim()[1] / 10) * 8
    right_upper.text(x, y, 
                     'FDR: %.2e\nEnrichment: %.2f'%(pval, fold),
                     fontsize = 10,
                    va = 'top', fontweight = 'bold')
    ## plot bars
    right_lower.add_patch(Rectangle((0, 0), N, 1, color = '#F8F8FF'))
    right_lower.bar(indices, 1, width = 50, color = dot_color)
    right_lower.spines['right'].set_visible(False)
    right_lower.spines['top'].set_visible(False)
#     right_lower.spines['bottom'].set_visible(False)
#     right_lower.spines['left'].set_visible(False)
    right_lower.set_ylabel('')
    right_lower.set_yticks([])
    right_lower.set_xticks([])
    right_lower.set_xlabel('Weight\n(Scaled Gene Expression & Correlation)',
                          fontsize = 12)
    right_lower.set_xlim(-100, N)
    fig.suptitle(title)
    plt.tight_layout()
    pdf.savefig(fig) if pdf else None
    if show_plot:
        plt.show()
    plt.close()
    if return_fig:
        return(fig)

def _cellpair_ES_plot_(geneSet, 
                  mHG_obj,
                  expList,
                  gxtList,
                  sender = '',
                  receiver = '',
                  figsize = (8, 3.5),
                  dot_color = '#1874CD',
                  curve_color = 'black',
                  title = '',
                  pdf = None,
                  show_plot = False,
                  return_fig = False
             ):
    ##clean
    plt.close() 
    ## geneset
    cgene = np.intersect1d(expList.index, gxtList.index)
    geneSet = np.intersect1d(geneSet, cgene)
    ## mHG obj
    N = mHG_obj.N
    indices = mHG_obj.indices
    pvals, folds = get_hypergeometric_stats(mHG_obj.N, mHG_obj.indices)
    try:
        fold, pval = mHG_obj.fold_enrichment, mHG_obj.pval
    except:
        fold, pval = 0, mHG_obj.pval
    ## figure
    fig = plt.figure(constrained_layout=True, figsize = figsize)
    subfigs = fig.add_gridspec(2, 2, width_ratios=[1,1],
                               height_ratios = [1,0.1],
                               hspace = 0, wspace = 0.05)
    left = fig.add_subplot(subfigs[:, 0])
    right_upper = fig.add_subplot(subfigs[0, 1])
    right_lower = fig.add_subplot(subfigs[1, 1], sharex = right_upper)
    ## scatter
    ## left scatter
    mx = np.median(expList)
    left.scatter(x = expList[cgene], y = gxtList[cgene], s = 1, 
               edgecolor = 'none', facecolor = 'lightgrey')
    ymin, ymax = left.axes.get_ylim()
    left.vlines(mx, ymin, ymax, color = 'black',
              linestyle = '--', linewidth = .7)
    xmin, xmax = left.axes.get_xlim()
    left.hlines(0, xmin, xmax, color = 'black',
              linestyle = '--', linewidth = .7)
    left.scatter(x = expList[geneSet], y = gxtList[geneSet], 
               color = dot_color, s = 10, alpha = .8)
    ## ann
    q1, q2, q3, q4 = _text_quadrant(expList, gxtList, geneSet)
    
    extend = 0
    left.text(x = xmax, y = ymax - extend,
            s = '{}/{}'.format(q1,len(geneSet)),
            color = dot_color, ha = 'right')
    left.text(x = xmin, y = ymax - extend,
            s = '{}/{}'.format(q2,len(geneSet)),
            color = dot_color, ha = 'left')
    left.text(x = xmin, y = ymin + extend,
            s = '{}/{}'.format(q3,len(geneSet)),
            color = dot_color, ha = 'left')
    left.text(x = xmax, y = ymin + extend,
            s = '{}/{}'.format(q4,len(geneSet)),
            color = dot_color, ha = 'right')
    if receiver != '':
        left.set_xlabel('Scaled Expression in %s'%receiver, fontsize = 12, labelpad = 10)
    else:
        left.set_xlabel('Scaled Expression', fontsize = 12, labelpad = 10)
    if sender != '' and receiver != '':
        left.set_ylabel('Aggregated Weight from Communications\nbetween %s and %s'%(sender, receiver), fontsize = 12)
    else:
        left.set_ylabel('Aggregated Weight', fontsize = 12)
#     sns.despine(trim = True, offset = 2, ax = left)
    left.spines['right'].set_visible(False)
    left.spines['top'].set_visible(False)
    ## enrich curve
    right_upper.plot(range(N+1), -np.log10(pvals), #marker = 1,
               color = curve_color, markersize = 5, linewidth = .7)
    right_upper.scatter(indices, -np.log10(pvals)[indices],
                        color = dot_color, s = 5, zorder = 3,
                      label = 'Pathway gene')

    right_upper.set_ylabel('Enrichment score', fontsize = 12)
    right_upper.spines['right'].set_visible(False)
    right_upper.spines['top'].set_visible(False)
#     right_upper.spines['bottom'].set_visible(False)
    right_upper.set_xticks([])
    ymin, ymax = right_upper.get_ylim()
    right_upper.text(0, ymin-(ymax*0.065), 'Highly associated', color = 'red')
    right_upper.text(N, ymin-(ymax*0.065), 'Lowly associated', 
               color = 'blue', ha = 'right')
    leg = right_upper.legend(fontsize = 10, frameon = False)
    right_upper.set_xlim(-100, N)
    ## label scores
    x = (N / 3) * 1.7
    y = (right_upper.get_ylim()[1] / 10) * 8
    right_upper.text(x, y, 
                     'FDR: %.2e\nEnrichment: %.2f'%(pval, fold),
                     fontsize = 10,
                    va = 'top', fontweight = 'bold')
    ## plot bars
    right_lower.add_patch(Rectangle((0, 0), N, 1, color = '#F8F8FF'))
    right_lower.bar(indices, 1, width = 50, color = dot_color)
    right_lower.spines['right'].set_visible(False)
    right_lower.spines['top'].set_visible(False)
#     right_lower.spines['bottom'].set_visible(False)
#     right_lower.spines['left'].set_visible(False)
    right_lower.set_ylabel('')
    right_lower.set_yticks([])
    right_lower.set_xticks([])
    right_lower.set_xlabel('Aggregated Weight\n(Scaled Gene Expression & Correlation)',
                          fontsize = 12)
    right_lower.set_xlim(-100, N)
    fig.suptitle(title)
    plt.tight_layout()
    pdf.savefig(fig) if pdf else None
    if show_plot:
        plt.show()
    plt.close()
    if return_fig:
        return(fig)

