![LOGO](./images/mebocost_logo.png | width=10)
## Welcome to use MEBOCOST, a python-based package to predict metabolite-based cell-cell communications by single-cell RNA-seq data of tissue samples.

### Version control
<li>Current release: 1.0.0</li>

### Installation
* download miniconda enviroment
```{bash}
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh && bash Miniconda3-latest-MacOSX-x86_64.sh

conda create -n mebocost python=3.8

conda activate mebocost
```
* download MEBOCOST package from github
```{bash}
git clone xxxx

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
#### To check whether it has been installed sucessfully, users can run in python:
```{python}
>>from mebocost import mebocost
```
#### if the mebocost can be imported successfully, you can continue to do analyses by mebocost!

### Tutorial of MEBOCOST

