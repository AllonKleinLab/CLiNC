# Cell-Lineage-from-Normalized-Covariance


Cell-Lineage-from-Normalized-Covariance (CLiNC) is a method to reconstruct developmental hierarchies from clonal barcoding data. The method is described in REF. Briefly, the model underlying CLiNC assumes that all barcodes are deposited as a synchronous moment in differentiation and that differentiation events are not directly coupled to cell division (as in asymmetric division). 

#### Algorithm overview
The input data is a matrix of barcode counts in across cell types. In principle these counts should represent numbers of cells (as opposed to numbers of sequencing reads). The output is an inferred cell type hierarchy and a list of putative tree violations. The only parameter is the false-discovery rate for detection of conformal symmetry violations (default 5%). The CLiNCs pipeline includes the following steps:

1. Calculate normalized covariance between each pair of cell types
2. Use neighbor-joining to iterative form a cell type hierarchy
3. Identify statistically significant deviations from conformal symmetry
4. Use symmetry violations to infer putative differentiation pathways that violate the hierarchy


## Installation ##


#### Install Miniconda

First, you will need to download and install the [Miniconda3 environment from Continuum Analytics](https://conda.io/miniconda.html)

Get the *Python 3.6* version for your operating system.  If you are on a Mac or Linux machine, run the bash installer. 

On Linux,

```sh
# downloads the anaconda installer
curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o "$HOME/miniconda3_latest.sh"
# makes the installer an executable
chmod +x $HOME/miniconda3_latest.sh
# run the installer
$HOME/miniconda3_latest.sh -b -p $HOME/miniconda3
```

On Mac OS X,

```sh
# downloads the anaconda installer
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o "$HOME/miniconda3_latest.sh"
# makes the installer an executable
chmod +x $HOME/miniconda3_latest.sh
# run the installer
$HOME/miniconda3_latest.sh -b -p $HOME/miniconda3
```

On Windows 10, 

First install the [Debian WSL](https://www.microsoft.com/en-us/p/debian/9msvkqc78pk6?activetab=pivot:overviewtab) from the Windows store.  Then, open up the app, which should drop you into a Bash terminal.  Run,

```sh
sudo apt-get update
sudo apt-get install git build-essential curl libxrender-dev libsm6 libglib2.0-0
# downloads the anaconda installer
curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o "$HOME/miniconda3_latest.sh"
# makes the installer an executable
chmod +x $HOME/miniconda3_latest.sh
# run the installer
$HOME/miniconda3_latest.sh -b -p $HOME/miniconda3
```

Then add the appropriate lines to your `.bashrc` file (or `.bash_profile` on a Mac)

```sh
cat >> ~/.bashrc << END
PATH=\$HOME/miniconda3/bin:\$PATH
END
source $HOME/.bashrc
```

#### Create virtual environment

Then create a new virtual environment called `clinc` using `conda` and activate it,

```sh
conda create -n "clinc" python=3.6 -y
conda activate clinc
```

#### Installing dependencies

(Run each of the following lines separately)
```
conda activate clinc
conda install ipykernel jupyter numpy scipy matplotlib
conda install -c etetoolkit ete3
conda install -c conda-forge fastcluster
python -m ipykernel install --user --name clinc --display-name "Python (clinc)"
```

## Usage ##

Download and ```cd``` into the CLiNC repository:
```
git clone https://github.com/AllonKleinLab/Cell-Lineage-from-Normalized-Covariance.git
cd Cell-Lineage-from-Normalized-Covariance
```

Start jupyter
```jupyter notebook ```
 and click on ```CLiNC pipeline.ipynb```
