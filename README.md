
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
tsfm -h
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

1. Recreating the function logos for the human tRNA data. This directory contains all of the aligned tRNA sequences used for analysis in Kelly et al. 2020.
   1. First we need to change into the example data directory
   ```shell
   cd Kelly2020_data
   ```
   
   2. To create single site and basepair function logos for the human tRNA data using the NSB entropy estimator we can use this command:
   ```shell
   tsfm -e NSB -x 5 -c tRNA_L_skel_Leish.sites74.struct.cove --logo -p 10 HOMO/HOMO
   ```
   3. Lets break this command down so we can understand the options   
      a. The below `-e` option sets the entropy estimator to NSB. This can also be set to `Miller` to use the Miller-Madow estimator.
      ```shell
      -e NSB
      ```
      b. The `-x` option indicates the maximum sample size for the calculation of the exact entropy correction. Correction values can be feasibly calculated for sample sizes up to ~16.
      ```shell
      -x 5
      ```
      c. The `-c` option provides the path to the file containing the secondary structure annotation in cove format. This is required for basepair function logos and can be provided in `inferenal`, `cove`, or `text` format. These format as described in below in Appendix section. 
      ```shell
      -c tRNA_L_skel_Leish.sites74.struct.cove
      ```
      d. The below option indicates that we want graphic output of function logos. If we excluded this option, tSFM would only produce the text output file.
      ```shell
      --logo
      ```
      c. The `-p` option indicates the number of processor cores we want to utilize during the calculations. Below we have indicated we want to use 10 cores.
      ```shell
      -p 10
      ```
      e. Lastly, the below argument provides the path and prefix to the alignment files containing our tRNA sequences. For function/ID/KLD logos tSFM expects alignments partitioned into "taxa" by the prefix and functional class by the postfix. For example, in the `HOMO` directory, we have the alignement file `HOMO_A.aln`, which the prefix `HOMO` indicates the "taxa" and the postfix `_A` indicates the alignment file contains sequences for the Alanine functional class.  
      ```shell
      HOMO/HOMO
      ```
2. To avoid having to manually enter the command to calculate function logos for each of the parasite clades we can take a shortcut by using bash loops. The below command will loop through each folder and produce function logos for each clade.
```shell
for d in */; do tsfm -e NSB -x 5 -c tRNA_L_skel_Leish.sites74.struct.cove --logo -p 10 $d${d%/}; done
```
3. Recreating ID/KLD logos and data table for bubble plots for `HOMO` versus `MAJOR` 
```shell
tsfm -c tRNA_L_skel_Leish.sites74.struct.txt -e Miller -x 5 --idlogos --kldlogos --bt MAJOR/MAJOR HOMO/HOMO
```

# Statistical significance testing
# Recreating the supplemental figure from the tSFM publication
# Appendix 
## Structure file
