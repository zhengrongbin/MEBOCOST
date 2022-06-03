<img src="./images/mebocost_logo.png" width="200" height="180" style='margin-left: auto; margin-right: auto;display: block;'></img>

## Welcome to use MEBOCOST: <I>Me</I>ta<I>bo</I>lic Cell-Cell <I>Co</I>mmunication Modeling by <I>S</I>ingle Cell <I>T</I>ranscriptome

### What is MEBOCOST and hoe does it work?
<p>MEBOCOST is a Python-based computational tool for inferring metabolite, such as lipid, mediated cell-cell communication events using single-cell RNA-seq data. Briefly, in the first step, MEBOCOST imputes the relative abundance of metabolites based on the gene expression of metabolic reaction enzymes. The genes of enzymes were collected from Human Metabolome Database (HMDB). Next, MEBOCOST identifies cell-cell metabolite-sensor communications between cell groups, in which metabolite enzymes and sensors were highly expressed in sender and receiver cells, respectively. Furthermore, MEBOCOST can infer communication associated pathways in receiver cells, which we try to link cellular mechanism of receivers in response to the communication events.</p>

#### The Flowchart of MEBOCOST
<p>Panel A: workflow for predicting cell-cell metabolic communication events taking scRNA-seq data as input. Panel B for pathway association inference for the significant communication events.</p>
<img src='./images/FlowChart.png' style='margin-left: auto; margin-right: auto;display: block;'></img>

### Version control
<p>We keep updating MEBOCOST!!!</p>
<li>Current release: 1.0.1</li>
- add a "min_cell_number" parameter to enable the exclusion of cell groups that do not have enough cell numbers
- fixed bug for notebook visualization for flow plot
<hr>

### Installation
* download and install miniconda enviroment (Users can skip this step if a python-based environment has been well-established)
```{bash}
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh && bash Miniconda3-latest-MacOSX-x86_64.sh

conda create -n mebocost python=3.8

conda activate mebocost
```
* download MEBOCOST package from github
```{bash}
git clone https://github.com/zhengrongbin/MEBOCOST.git

cd MEBOCOST
```
* install requirements
```{bash}
pip install -r requirements.txt
```
* install MEBOCOST
```{bash}
python setup.py install
```

### To re-install and upgrade by pip if you have existing MEBOCOST installed
```
pip uninstall mebocost
pip install git+https://github.com/zhengrongbin/MEBOCOST.git --upgrade
```


#### To check whether it has been installed sucessfully, users can run in python:
```{python}
>>from mebocost import mebocost
```
#### if the mebocost can be imported successfully, you can continue to do analyses by mebocost!

### Tutorial for MEBOCOST

<li><a href='./Demo_Communication_Prediction.ipynb' target='_blank'>Prediction of cell-cell metabolic communication by scRNA-seq data</a></li>
<li><a href='./Demo_Pathway_Inference.ipynb' target='_blank'>Inference of cell-cell metabolic communication associated pathways in receiver cells</a></li>

### Cite us
<p>Please cite us at <a href='https://www.biorxiv.org/content/10.1101/2022.05.30.494067v1' target='_blank'>bioRxiv</a> if you find MEBOCOST is useful to your project.</p>

### Contact
Rongbin.Zheng@childrens.harvard.edu

or

Kaifu.Chen@childrens.harvard.edu

<hr>
Copy Right @ Kaifu Chen Lab @ Boston Childrens Hospital / Harvard Medical School
