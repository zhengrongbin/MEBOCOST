## Manual for MEBOCOST software
<p>We explain each parameters in MEBOCOST functions</p>

#### 1.1 Initial MEBOCOST project
```{python}
>> from mebocost import mebocost
>> mebo_obj = mebocost.create_obj(
                        adata = adata,
                        group_col = ['celltype'],
                        met_est = 'mebocost',
                        config_path = './mebocost.conf',
                        exp_mat=None,
                        cell_ann=None,
                        species='human',
                        met_pred=None,
                        met_enzyme=None,
                        met_sensor=None,
                        met_ann=None,
                        scFEA_ann=None,
                        compass_met_ann=None,
                        compass_rxn_ann=None,
                        gene_network=None,
                        gmt_path=None,
                        cutoff_exp=0,
                        cutoff_met=0,
                        cutoff_prop=0.25,
                        sensor_type=['Receptor', 'Transporter', 'Nuclear Receptor'],
                        thread=8
                        )
```
* Params for creat_obj
-------
| Parameter | Default | Data type | Description
| :------------------------ |:-------------:| :-------------|
| exp_mat | None | pandas data frame | single cell expression matrix, rows are genes, columns are cells, this is exclusive to 'adata' |
| adata | None | scanpy object (adata) | scanpy adata object, the expression will be extracted, 'adata' is exclusive to 'exp_mat' |
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