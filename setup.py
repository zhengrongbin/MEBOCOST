"""lisa: a bioinformatics software
epigenome analysis to rank TFs from gene set
"""
import os
import setuptools
import configparser

def _make_config(conf_path, workdir=os.getcwd()):
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
                path = config[firstLevel][secondLevel][:config[firstLevel][secondLevel].index('#')-1].rstrip()
                config[firstLevel][secondLevel] = os.path.join(workdir, path)
            else:
                path = config[firstLevel][secondLevel]
                config[firstLevel][secondLevel] = os.path.join(workdir, path)
    ## re-write
    cf_new = configparser.ConfigParser()
    for firstLevel in config.keys():
        cf_new.add_section(firstLevel)
        for secondLevel in config[firstLevel]:
            cf_new.set(firstLevel, secondLevel, config[firstLevel][secondLevel])
    with open('mebocost.conf', 'w') as f:
        cf.write(f)
    return(config)

## setup
def main():
  setuptools.setup(name="mebocost", 
                  version="1.0.0",
                  description="a python-based method to predict metabolite mediated cell-cell communication",
                  author='Rongbin Zheng, Kaifu Chen',
                  author_email='Rongbin.Zheng@childrens.harvard.edu',
                  url='http://liulab.dfci.harvard.edu/',
                  # scripts=glob('mebocost/*'),
                  zip_safe=True,
                  package_dir={"": "src"},
                  packages=setuptools.find_packages(where="src"),
                  classifiers=[
                      'Environment::Console',
                      'Operating System:: POSIX',
                      "Programming Language :: Python :: 3",
                      "Topic :: Scientific/Engineering :: Bio-Informatics"],
                  keywords='Metabolism',
                  license='OTHER'
  )
if __name__ == '__main__':
    ## change mebocost.conf to absolute path
    _make_config(conf_path = './src/mebocost.conf')
    ## setup
    main()
