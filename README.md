
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

Kelly, P., F. Hadi-Nezhad, D. Y. Liu, T. J. Lawrence, R. G. Linington, M. Ibba, and D. H. Ardell. 2020. Targeting tRNA-synthetase interactions towards novel therapeutic discovery against eukaryotic pathogens. PLOS Neglected Tropical Diseases 14: e0007983.
