
# Installation
1. Download tSFM and example data using the below bash command or the zip link: 
```bash
git clone https://github.com/tlawrence3/tSFM.git
```
2. Change into the `tSFM` directory
```bash
cd tSFM
```
3. We recommend using [anaconda](https://www.anaconda.com/products/individual) and conda enviroments for installing and using tSFM. The following commands will create a conda environment with all the required packages and activate it.
```bash
conda create -n tSFM python=3.8 pandas numpy cython pytest-runner pytest scipy pasty mpmath statsmodels
conda activate tSFM
```
4. Now we can install tSFM and run a basic test the installation with the below commands:
```bash
python setup.py install
tSFM -h
```

# Running more extensive tests
If you want to run more extensive testing that we utilize during development you can use the below commands:
```shell
pip install pytest
pip install -r requirements.txt
python setup.py build_ext --inplace
python -m pytest tests/
```

# Quickstart tutorial
As a quick introduction to the functionality of of tSFM we will be utilizing the data from: 

[Kelly, P., F. Hadi-Nezhad, D. Y. Liu, T. J. Lawrence, R. G. Linington, M. Ibba, and D. H. Ardell. 2020. Targeting tRNA-synthetase interactions towards novel therapeutic discovery against eukaryotic pathogens. PLOS Neglected Tropical Diseases 14: e0007983.](https://doi.org/10.1371/journal.pntd.0007983)

1. Recreating the function logos for the human tRNA data.
   1. First we need to change into the example data directory
   ```shell
   cd Kelly2020_data
   ```
   This directory contains all of the aligned tRNA sequences used for analysis in Kelly et al. 2020.
   1. To create single site and basepair function logos for the human tRNA data using the NSB entropy estimator we can use this command:
   ```shell
   tsfm -e NSB -x 5 -c tRNA_L_skel_Leish.sites74.struct.cove --logo HOMO/HOMO
   ```
# Recreating the supplemental figure from the tSFM publication
