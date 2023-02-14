## Directory

    .
    ├── data  
    |    ├── tpm_sample_genes.csv
    |    ├── mm10 / genes / genes.gtf
    |    └── loom
    |         ├── E11A_all.loom  
    |         ├── E11Y_all.loom  
    |         ├── E12F_all.loom  
    |         └── E14F_all.loom  
    ├── Fig5A                
    |    └── fig5a.ipynb
    ├── Fig5B                
    |    └── fig5b.ipynb
    ├── TableS1    
    |    ├── deseq_right_vs_left.R
    |    ├── deseq_right_end_vs_left_end.R
    |    ├── x_coordinate_right_vs_left.csv
    |    └── x_coordinate_right_end_vs_left_end.csv
    ├── poetry.lock  
    ├── pyproject.toml  
    └── README.md  



## Setup
### pyenv install (if needed)
```sh
$ git clone git://github.com/yyuu/pyenv.git ~/.pyenv  
$ echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bash_profile  
$ echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bash_profile  
$ echo 'eval "$(pyenv init --path)"' >> ~/.bash_profile  
$ source ~/.bash_profile  
```

### pyenv setup  
```sh
$ pyenv install 3.7.0  
$ cd ~/hsc_independent #project directory name  
$ pyenv local 3.7.0  
```


### poetry install
```sh
$ pip3 install poetry==1.1.8  
$ poetry config virtualenvs.in-project true
```


### package install
```sh
$ cd ~/hsc_independent #project directory name 
$ poetry install  
$ poetry run jupyter lab  
```  


### jupyter remote access (if needed)
https://jupyter-notebook.readthedocs.io/en/stable/public_server.html


### Failing to install python-igraph (if needed)  
Setup your computer environment for package compiling

(Ubuntu)
```sh
$ sudo apt update gcc build-essential python3-dev libxslt-dev libffi-dev libssl-dev libxml2 libxml2-dev zlib1g-dev
```