
# Installation
1. Download tSFM and example data using the below bash command or the zip link: 
```bash
git clone https://github.com/tlawrence3/tSFM.git
```
2. Change into the `tSFM` directory
```bash
cd tSFM
```
3. We recommend using
[anaconda](https://www.anaconda.com/products/individual) and conda
enviroments for installing and using tSFM. The following commands will
create a conda environment with all the required packages and activate
it. The last two commands install scikit-gof from its github repo
including two unmerged pull requests to fix installation errors and a
runtime warning. scikit-gof implements the Anderson-Darling Test used
in the Peaks-over-Threshold-GPD based p-value estimator.
```bash
conda create -n tSFM python=3.8 pandas numpy cython pytest-runner pytest scipy patsy mpmath statsmodels
conda activate tSFM
pip install git+https://github.com/wrwrwr/scikit-gof.git@952674f186c70077deb703b3ea39ee9bdd58a0aa
pip install git+https://github.com/wrwrwr/scikit-gof.git@b441372ecdbe110313720eea45d7e4833a279b24
```
4. Now we can install tSFM and run a basic test the installation with the below commands:
```bash
python setup.py install
tsfm -h
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
tSFM implements statisitical significance testing using permutation based null distributions and corrects multiple testing using FDR and FWER methods. The `-P` option indicates the number of permutations generated for building the null distributions and the `-C` option indicates the method for multiple testing correction.    

1. To calculate the significance of CIFs stack-heights with 100 permutations for the humans tRNA data using the NSB entropy estimator we can use this command:
```shell
tsfm -e NSB -x 5 -c tRNA_L_skel_Leish.sites74.struct.cove --logo -T stacks -P 100 HOMO/HOMO
```
   The `-T` option will set the significance test for only CIF stack-heights. This can also be set to `letters` to calculate the significance for only CIF letter- heights. the default is for both.
   ```shell
   -T stacks 
   ```
   The `-P` option will set the number of permutations for significance test of CIFs to 100. 
   ```shell
   -P 100
   ```

2. To calculate the significance of KLD-values with 100 permutations for clade humans and clade MAJOR against each other we can use this command:
```shell
tsfm --kldperms 100 -c tRNA_L_skel_Leish.sites74.struct.cove -e NSB -x 5 HOMO/HOMO MAJOR/MAJOR
```
   The `--kldperms` option will set the number of permutations to compute significance of Kullback-Leibler Divergences (KLD values) to 100
   ```shell
   --kldperms 100
   ```

3. To calculate the significance of ID-values for clade HOMO against clades MAJOR and ENRIETTII with 100 permutations we can use this command:
```shell
tsfm --idperms 100 -c tRNA_L_skel_Leish.sites74.struct.cove -e NSB -x 5 --clade HOMO HOMO/HOMO MAJOR/MAJOR ENRIETTII/ENRIETTII
```
   The `--idperms` will set the number of permutations to compute significance of Information Differences (ID values) to 100
   ```shell
   --idperms 100
   ```
   The `--clade` option will contrast two clades MAJOR and ENRIETTII against HOMO
   ```shell
   --clade HOMO
   ```

# Recreating the supplemental figure from the tSFM publication
Warning: this will take ~14 hours on 24 cores
```shell
tsfm --kldperms 10000 --idperms 10000 -c tRNA_L_skel_Leish.sites74.struct.cove -e NSB -x 5 --clade HOMO HOMO/HOMO ENRIETTII/ENRIETTII MAJOR/MAJOR
```

# Usage
```
file_prefix             One or more paths/file-prefix strings corresponding to
                        sets of input files compiled for a single clade in
                        clustalW format. Input files should be named
                        <path>/<prefix>_<functional-class>.<extension>, where
                        <functional_class> is a single letter

optional arguments:
  -h, --help            show this help message and exit
  -i INFERNAL, --infernal INFERNAL
                        Use secondary structure file INFERNAL, required to
                        calculate functional information of base-pair
                        features, in Infernal format
  -c COVE, --cove COVE  Use secondary structure file COVE, required to
                        calculate functional information of base-pair
                        features, in COVE format. Example: 
                        #=CS >>>>>>>..>>>>...........<<<<.>>>>>.......<<<<<.....>>>>>....
                        #=CS ...<<<<<<<<<<<<.
  -t TEXT, --text TEXT  Use secondary structure file TEXT, required to
                        calculate functional information of base-pair
                        features, is in text format. Example:
                        "A:0,72,1,71,2,70,3,69,4,68,5,67,6,66
                        D:9,25,10,24,11,23,12,22
                        C:27,43,28,42,29,41,30,40,31,39
                        T:49,65,50,64,51,63,52,62,53,61"
  -s, --single          Do not calculate functional information for paired
                        features. Calculate for single-site features only.
  -n, --nosingle        Do not calculate functional information for single-
                        site features. Calculate for paired-site features
                        only.
  -V, --version         show program's version number and exit
  -p PROCESSES, --processes PROCESSES
                        Set the maximum number of concurrent processes.
                        Default is the number of cores reported by the
                        operating system.
  -e {NSB,MM}, --entropy {NSB,MM}
                        Use entropy estimator ENTROPY when conditional sample
                        sizes exceed maximum for exact calculation. If value
                        is "NSB", use Nemenman-Shafee-Bialek estimator. If
                        value is "MM", use Miller-Madow estimator. Default is
                        NSB
  -x EXACT, --exact EXACT
                        Maximum conditional sample size to exactly calculate
                        entropy correction. Default = 5
  -l, --logos           Visualize feature information using function/inverse
                        function logos in extended postscript format.
  -v, --inverse         Additionally calculate functional information for
                        inverse features/logos (features under-represented in
                        specific functional classes)
  -P PERMUTATIONS, --permutations PERMUTATIONS
                        Calculate the significance of CIFs by a permutation
                        test, with a number of permutations equal to
                        PERMUTATIONS (an integer). Default is to not calculate
                        significance of CIFs.
  -C, --correction {bonferroni,sidak,holm,holm-sidak,simes-hochberg,hommel,BH,BY,GBS}
                        Specify a method for multiple test correction for
                        significance calculations: bonferroni, sidak, holm,
                        holm-sidak, simes-hochberg, hommel, BH (Benjamini-
                        Hochberg FDR), BY (Benjamini-Yekutieli FDR) or GBS
                        (Gavrilov-Benjamini-Sarkar FDR). Default is BH
  -T {stacks,letters,both}
                        Test the significance of only CIF stack-heights, only
                        CIF letter-heights, or both. Default is both.
  -I, --idlogos         Compute Information Differece logos for each pair of
                        clades
  -K, --kldlogos        Compute Kullback-Liebler Divergence logos for each
                        pair of clades
  -B, --bubbles         Compute input table for structural bubble-plots (to be
                        computed in R) like those appearing in Kelly et al.
                        (2020).
  --kldperms KLDPERMS   Set the number of permutations to compute significance
                        of Kullback-Leibler Divergences. Default is to not
                        calculate signifiance.
  --idperms IDPERMS     Set the number of permutations to compute significance
                        of Information Differences (SLOW). Default is to not
                        calculate signifiance.
  -J, --JSD             Produce pairwise distance matrices between function
                        logos for different taxa based on Jensen-Shannon
                        Divergences
  --clade CLADE         Contrast clade CLADE against all others. CLADE should
                        be one of the file-prefix strings passed as a required
                        argument to the program, stripped of its path. Default
                        is to compute contrasts for all pairs of clades.
```
