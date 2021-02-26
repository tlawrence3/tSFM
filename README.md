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
it. 
```bash
conda create -n tSFM python=3.8 pandas numpy cython pytest-runner pytest scipy patsy mpmath statsmodels
conda activate tSFM
```
4. Now we can install tSFM and run a basic test the installation with the below commands:
```bash
python setup.py install
tsfm -h
```

# Quickstart Tutorial
tSFM may be used to analyze any protein or RNA families. As a quick introduction to the functionality of tSFM we will be utilizing data from: 

[Kelly, P., F. Hadi-Nezhad, D. Y. Liu, T. J. Lawrence, R. G. Linington, M. Ibba, and D. H. Ardell. 2020. Targeting tRNA-synthetase interactions towards novel therapeutic discovery against eukaryotic pathogens. PLOS Neglected Tropical Diseases 14: e0007983.](https://doi.org/10.1371/journal.pntd.0007983)


Let us recreate function logos for the human tRNA data provided by this study.

1. First we need to change into the example data directory. This directory contains all of the aligned tRNA sequences used for analysis in Kelly *et al.* 2020.
    ```shell
    cd Kelly2020_data
    ```
   
2. To create single site and basepair RNA function logos for the human tRNA data we can use the command:
    ```shell
    tsfm -e NSB -x 1 -c tRNA_L_skel_Leish.sites74.struct.cove --logo -p 10 HOMO/HOMO 
    ```
   
    Let us break this command down so we can understand the options:
   
      a. The below `-e` option sets the entropy estimator to the Bayesian estimator NSB (which is the default). This can also be set to the older and more inaccurate Miller-Madow estimator with argument `Miller`.

      ```shell
	 -e NSB
      ```
	
      b. The `-x` option indicates the maximum sample size for the calculation of the exact entropy correction (larger sample sizes use the estimator set with the `-e` option). Correction values can feasibly be calculated for sample sizes up to ~16, but the NSB estimator is preferred except for the smallest sample sizes.

      ```shell
	   -x 1
      ```
	
      c. The `-c` option provides the path to the file containing the
      secondary structure annotation in cove format for RNA input
      files. This is required for basepair function logos and can be
      provided in `infernal`, `cove`, or `text` formats. These formats
      are described in below in the Appendix section. The secondary
      structure annotation must conform to and encode the RNA
      structural alignment common to the input alignment files and
      usually must be prepared together with them.

      ```shell
        -c tRNA_L_skel_Leish.sites74.struct.cove
      ```
	
      d. The below option indicates that we want graphical output of function logos in EPS format. If we excluded this option, tSFM would only produce the text output file with statistics on tRNA CIFs.

      ```shell
        --logo
      ```
	
      e. The `-p` option indicates the number of processor cores we want to utilize during the calculations. Below we have indicated we want to use 10 cores.

      ```shell
        -p 10
      ```
	
      f. Lastly, the below argument provides the path and clade prefix
      to the input alignment files containing our tRNA
      sequences. Multiple sets of input files for multiple clades can
      be processed at the same time by providing multiple
      arguments. tSFM expects input alignment files in ClustalW
      format, all named with a common prefix corresponding to their
      clade followed by a distinct one-letter symbol for the
      functional class of sequences in the alignment, separated by an
      underscore character.  For example, in the `HOMO` directory, we
      have multiple alignment files including `HOMO_A.aln`, in which
      the prefix `HOMO` indicates the clade and the postfix `_A`
      indicates the alignment file containing sequences for the
      tRNA-Ala functional class.

      ```shell
        HOMO/HOMO
      ```
	
3. For each clade on input, tSFM produces a text output file in the
   same directory where the program is run. This command produces the
   output file `HOMO_CIFs.txt`, more generally
   `<clade-prefix>_CIFs.txt`. The text output file has two sections,
   one for paired-site features and one for single-site features.
   Examples of records from each of these two sections is shown below.

   #bp |	coord	| state	| N |	info	|p-value|	BH|class:height:p-value:BH
   :-: | :-: | :-: | :--: | :-: | :-: | :-: | :-: 
   bp:|	(0, 72)|	AU|	14|	2.789|	NA|	NA|	X:0.861:NA:NA  L:0.108:NA:NA  K:0.032:NA:NA 

   The first record in the output file is for feature `(0, 72)AU`: the
   number of sequences sharing this features is in column `N`, the
   total information (the total height of the stack (in bits) from
   function logo of feature `AU` at coordinate `(0,72)`) is in column
   `info`, the p-value of the calculated info is in column `p-value`
   and the multiple-test corrected significance of the p-value is in
   the next column (this column is named with the method used for
   multiple test correction for significance calculations). Also, the
   last column shows a list of information for each tRNA class at
   feature `(0, 72)AU` with the format
   `<class>:<symbol-height>:<pvalue>:<BH>` which refers to the tRNA
   functional class single letter symbol, the height of the symbol
   within the stack, the p-value and the corrected significance value,
   respectively.

   A record for the single site is similar to the paired-site except that the first column is `ss` instead of `bp`.

   #ss |	coord	| state	| N |	info	|p-value|	BH|class:height:p-value:BH
   :-: | :-: | :-: | :--: | :-: | :-: | :-: | :-: 
   ss: |	0	| C	| 13|	4.131	|NA	|NA	|Y:1.000:NA:NA


   The values for the p-values and the significance are NA because we did not calculate p-values in the command described above. See the section `Statistical significance testing for CIFs` below, which describes options for p-value calculations.
         
   
4. To avoid having to manually enter the names and prefixes for each
   of the parasite clades we can take a shortcut by using bash
   loops. The below command will loop through each folder and produce
   function logos for each clade in the running directory.

   ```shell
   for d in */; do tsfm -e NSB -x 5 -c tRNA_L_skel_Leish.sites74.struct.cove --logo -p 10 $d${d%/}; done
   ```

5. To create ID/KLD logos and the data table for the bubble plots in
   Kelly *et al.* 2020. for clades `HOMO` versus `MAJOR`, we can run
   the command below.  Note that option --bubbles or -B requires
   designation of a specific clade to use as a reference clade (for
   gains and losses), using option --clade.

   ```shell
   tsfm -c tRNA_L_skel_Leish.sites74.struct.cove -e MM -x 5 --idlogos --kldlogos -B --clade HOMO MAJOR/MAJOR HOMO/HOMO
   ```

6. The text output of this command will be two files
   `F_HOMO_B_MAJOR_Table.txt` and `F_MAJOR_B_HOMO_Table.txt` which can
   be used to create the bubble plots. Here "F" stands for
   "Foreground" and "B" for "Background." An example of a record from
   the output text file is shown below.
        
   |aa | coord | state | fbits | fht | gainbits | gainfht | lossbits | lossfht | convbits | convfht | x | y | sprinzl |
   |:-: | :-: | :-: | :--: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
   |K | 1 | A | 2.7602 | 0.0 |	0.7285 |	0.0 |	0.0	| 0.0	| 0.2913|	0.0	|0.0|	0.0|	0.0 |

   This example shows at feature 1A of ID logo: the total height of
   stack-bar in column `gainbits`, and the height of symbol K in this
   stack-bar in column `gainfht`. Also, it shows at feature 1A of KLD
   logo: the total height of stack-bar in column `convbits` and the
   height of symbol K in this stack-bar in column `convfht`. The
   columns `x`, `y` and `sprinzl` are set to 0 and will be filled
   later prior to creating the bubble plots by mapping each feature to
   the tRNA sprinzl coordinates.

# Statistical Significance Testing for CIFs
tSFM implements statisitical significance testing for CIFs and function logos using permutation-based null distributions and corrects for multiple tests using multiple user-determined FDR or FWER methods. The `-P` option indicates the number of permutations generated for building the null distributions and the `-C` option indicates the method for multiple testing correction.    

1. To calculate the significance of CIFs stack-heights with 100 permutations for the humans tRNA data using the NSB entropy estimator we can use this command:
```shell
tsfm -e NSB -x 5 -c tRNA_L_skel_Leish.sites74.struct.cove --logo -T stacks -P 100 HOMO/HOMO
```
   The `-T` option will set the significance test for only CIF total stack-heights. This can also be set to `letters` to calculate the significance for only CIF letter-heights within stacks. The default is to calculate significance for both stacks and letters, but for many applications and comparisons, we recommend to test only `stacks`. The difference plays into the number of comparisons or tests requiring correction, where generally fewer is better. If one is interested in only which CIFs are significant, use `stacks.` If one wants to test specific functional associations, use `letters,` otherwise use `both`.

   Similarly, you can restrict CIF calculation and significance
   testing to only single-site features with the `--single` option or
   to only paired-site features with the `--nosingle` option. The `-P`
   option will set the number of permutations for significance test of
   CIFs to 100.

# Statistical Significance Testing for CIF Divergences 
tSFM implements statistical significance testing for CIF divergences, specifically Information Difference (ID) logos and Kullback-Leibler Divergence (KLD) logos.
using permutation-based null distributions. P-values can be calculated
by one of three algorithms called `ECDF`, `ECDF_pseudo`, and `GPD` chosen by the user through
the option `-m`.

The algorithm `ECDF_pseudo` is the Monte Carlo permutation-based
approximation calculated according to the formula (number of
exceedances + 1)/(number of permutations + 1), where `exceedances`
means the fraction of permutation replicates yielding statistics as
large or larger than the observed divergence. Algorithm `ECDF` is
similar to `ECDF_pseudo` except that when the number of permutation
exceedances reaches the value of option `--exceedances` (with default
10) the p-value will be calculated according to the formula (number of
exceedances)/(number of permutations) and returned, otherwise when
there are fewer exceedances, the `ECDF_pseudo` formula is
used. Algorithm `GPD` is implemented according to the
"Peaks-over-Threshhold" (PoT) tail-approximation approach described as
algorithm APPROXIMATE in the manuscript, which includes the `ECDF` and
`ECDF_pseudo` algorithms as special cases. The default method is
`GPD`. tSFM also calculates the confidence intervals for all the
p-values except for the ones with zero exceedances calculated with
pseudocounts.

1. To calculate the significance of KLD divergence values with 100 permutations for the contrast of the human and MAJOR clades, we can choose among the three algorithms of GPD, ECDF and ECDF_pseudo as in the three examples below:

    a. Method `ECDF_pseudo` 
    ```shell
    tsfm --kldperms 100 -m ECDF_pseudo -c tRNA_L_skel_Leish.sites74.struct.cove HOMO/HOMO MAJOR/MAJOR
    ```
    The `--kldperms` option will set the number of permutations to compute significance of KLD values to 100.
    
    
    b. Method `ECDF`.
    ```shell
    tsfm --kldperms 100 --exceedances 5 -m ECDF --alpha 0.03 -c tRNA_L_skel_Leish.sites74.struct.cove HOMO/HOMO MAJOR/MAJOR
    ```
    The `--alpha` option will set the significance level to compute the confidence interval of pvalues. Its default is 0.05.
    
    
    c. Method `GPD`. In addition to the previous options, there are options `--targetperms` and `--peaks` specific to method `GPD`. Options `targetperms` and `peaks` are referred to as variables T and U, respectively in the algorithm APPROXIMATE. Also the option `exceedances` is referred to as parameter S and is used for both methods `ECDF` and `GPD`. 
    ```shell
    tsfm --kldperms 100 -m GPD --targetperms 70 --exceedances 5 --peaks 50 -c tRNA_L_skel_Leish.sites74.struct.cove HOMO/HOMO MAJOR/MAJOR
    ```
     i. The default value for option `targetperms` is 500. The value of the option `targetperms` should be less than the maximum permutation number indicated with option `--kldperms` or `--idperms`.
       
     ii. The default value of `exceedances` is 10. This number also needs to be less than the maximum permutation number and need not be (much) larger than 10, which is a standard rule-of-thumb for estimation of binomial proportions.
       
     iii. The default for option `peaks` is 250; However in the algorithm APPROXIMATE the peak will be set to the minimum of 250 and one-third of the permutations. The value of option peaks needs to be less than the maximum permutation number.
         
2. The output of KLD and ID logo significance from the examples described above will be two text files named `KLD_HOMO_MAJOR_stats.txt` and `KLD_MAJOR_HOMO_stats.txt`. 
   An example of a record from the output text file is shown below. This record shows the significance of the KLD statistic at feature U2 along with other information at this feature including: confidence interval (with level determined by option `--alpha`) in columns `CI.Lower` and `CI.Upper`, multiple-test adjusted p-value in column `Adjusted-P`, number of permutations with which the p-value is calculated in column `Permutations`, the method used for calculating the p-value in column `P-Val-Method` which can take the values: `p_ecdf`,  `p_ecdf_with_pseudo`, `p_ecdf_with_pseudo (p_gpd=0)` and `p_gpd`. If the p-value is calculated with GPD, the parameters of the GPD calculation will be shown in columns `GPD-shape`, `GPD-scale` and `Peaks` describing the maximum likelihood estimated parameters of the GPD distribution and the number of peaks over threshold. Also the column `ADtest-P-val` shows the pvalue of the goodness-of-fit test of the extreme permutation values to the GPD distribution.      
    
    Coord|State|Statistic|Sample-Sz-Back|Sample-Sz-Fore|P-value|CI.Lower|CI.Upper|Adjusted-P|Permutations|P-Val-Method|GPD-shape|GPD-scale|Peaks|ADtest-P-val|Freqs-Back|Freqs-Fore
    :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |:-:| :-: | :-: |:-:
    2|U|0.75|30|72|2.46e-14|0.0|3.22e-05|7.07e-13|70.0|p_gpd|0.15|0.00037|13.0|0.97|D13E16Y1|D24E24N24s

    
3. To calculate the significance of ID-values of single-site features only, for clade HOMO against clades MAJOR and ENRIETTII with 100 permutations using the method GPD we can use this command (note that running this command can take about 30 minutes.)
```shell
tsfm --idperms 100 -m GPD --targetperms 50 --exceedances 5 --peaks 50 --single -e NSB -x 5 --clade HOMO HOMO/HOMO MAJOR/MAJOR ENRIETTII/ENRIETTII
```
    a) The option `--idperms` will set the number of permutations to compute significance of ID values to 100
    ```shell
    --idperms 100
    ```
   
    b) The option `--clade` will contrast two clades MAJOR and ENRIETTII against HOMO
    ```shell
    --clade HOMO
    ```
   
    c) The option `--targetperms` will set the permutation number at which we attempt to fit the extreme permutation values to GPD.
    ```shell
    --targetperms 50
    ```

4. To calculate the significance of KLD-values of paired-site features only, for clade HOMO against clades MAJOR with 100 permutations using the method GPD use the command below. The consensus secondary structure of tRNAs indicated by option -c is required  to calculate functional information of base-pair features. 
```shell
tsfm --kldperms 100 -m GPD --targetperms 50 --exceedances 5 --peaks 50 --nosingle -c tRNA_L_skel_Leish.sites74.struct.cove --clade HOMO HOMO/HOMO MAJOR/MAJOR
```

# Usage and Options to tSFM
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
  -m {GPD,ECDF_pseudo,ECDF}, --pmethod {GPD,ECDF_pseudo,ECDF}
                        Set the p-value calculation algorithm for KLD/ID CIF Divergence significance calculations.
                        If value is "GPD", p-values for large divergences will be estimated by the Peaks-over-
                        Threshold method based on the Generalized Pareto Distribution, those for small divergences
                        with E exceedances (default 10) will be calculated as a binomial proportion (ECDF method)
                        and p-values for divergences that can't be estimated by either of those methods will be
                        estimated as a binomial proportion with pseudo-counts (ECDF_pseudo method). If value is
                        "ECDF", all p-values will be calculated by ECDF method or ECDF_pseudo method. If value is
                        "ECDF_pseudo", all the p-values will be calculated using ECDF_pseudo method. Default is GPD
  --targetperms TARGETPERMS
                        Set the initial target number of permutations at which GPD-pvalue will be initially
                        calculated. Default is 500.
  --exceedances EXCEEDANCES
                        Set the number of exceedances for which ECDF-based p-values will be calculated without
                        pseudo-counts. Default is 10.
  --peaks PEAKS         Set the number of Peaks-over-Threshold to use to initially estimate GPD. The actual value
                        used will be the minimum of this value and one-third of permutations. Default is 250.
  --alpha ALPHA         Set the significance level to compute the confidence interval of pvalues. Default is 0.05

```

# Recreating the Figure and Supplemental Figures from the tSFM Publication
1. Creating the supplemental figure <KLD-Significance_ENRIETTII_MAJOR.eps>
    
    1. Create the table of statistics for KLD logos of clade `ENRIETTII` vs `HOMO` using the pvalue-calculation-method `GPD` with maximum permutation number = 10000, target permutation number 500, and calculate the Confidence Interval of p-values with 5% significance level. The running time on an Intel Core i7 Dell XPS was ~34 minutes. The outputs are two text files: KLD_ENRIETTII_HOMO_stats.txt and KLD_HOMO_ENRIETTII_stats.txt 
    ```shell
    mkdir GPD
    cd Kelly2020_data
    tsfm --kldperms 10000 -m GPD -c tRNA_L_skel_Leish.sites74.struct.cove HOMO/HOMO ENRIETTII/ENRIETTII
    mv KLD_ENRIETTII_HOMO_stats.txt KLD_HOMO_ENRIETTII_stats.txt GPD
    ```
    
    2. Create the table of statistics for KLD logos of clade `MAJOR` vs `HOMO` using the pvalue-calculation-method `GPD` with maximum permutation number = 10000, target permutation number 500, and calculate the Confidence Interval of p-values with 5% significance level. The running time on an Intel Core i7 Dell XPS was ~52.6 minutes. The outputs are two text files: KLD_MAJOR_HOMO_stats.txt and KLD_HOMO_MAJOR_stats.txt
    ```shell
    tsfm --kldperms 10000 -m GPD -c tRNA_L_skel_Leish.sites74.struct.cove HOMO/HOMO MAJOR/MAJOR
    mv KLD_MAJOR_HOMO_stats.txt KLD_HOMO_MAJOR_stats.txt GPD
    ```
    
    3. Create the figure using the KLD tables in GPD folder by running this R script:
    ```R
    library(latex2exp)
    library(ggplot2)
    library(ggpubr)
    #setwd("PATH/TO/WHERE/THE/FOLDERS/GPD/ECDF/ECDF_pseudo/ARE/LOCATED") # make sure to set your working directory to where the input data folders GPD, ECDF and ECDF_pseudo are located.
    spetralc <-c("#D53E4F" ,"#F46D43" ,"#FDAE61" , "#66C2A5" ,"#3288BD","#313695")
    breaksarray <- c(1,2,4,10,30,80,200,500)
    legendbreaks <- logb(breaksarray, base = 2)
    legendlables <- as.character(breaksarray)
    KLD_paths <-
      c(
        "GPD/KLD_MAJOR_HOMO_stats.txt",
        "GPD/KLD_HOMO_MAJOR_stats.txt",
        "GPD/KLD_ENRIETTII_HOMO_stats.txt",
        "GPD/KLD_HOMO_ENRIETTII_stats.txt"
      )
    xlabels <- c("Kullback-Leibler Divergence","Kullback-Leibler Divergence","Kullback-Leibler Divergence","Kullback-Leibler Divergence")
    clades <- c("MAJOR (n=664 sequences) against Homo (m=431 sequences)", "HOMO against MAJOR","ENRIETTII (n=160 sequences) against Homo     (m=431 sequences)", "HOMO against ENRIETTII")
    plotls <- list()
    for (i in 1:length(KLD_paths)) {
      signTable <- read.table(KLD_paths[i], header = TRUE, sep = "\t")
      signTable <- signTable[signTable$Sample.Sz.Back != 0 & signTable$Sample.Sz.Fore != 0 & signTable$Statistic != 0, ]
      signTable$harmonicmean <- 0
      for (j in 1:nrow(signTable)) {
        signTable$harmonicmean[j] = 1 / mean(1 / c(signTable$Sample.Sz.Back[j], signTable$Sample.Sz.Fore[j]))
      }
      p <- ggplot(signTable, aes_string( x = round(signTable$Statistic, 4), y = -log(signTable$P.value, base = 2),colour =         logb(signTable$harmonicmean, 2))) + 
        geom_point(aes_string( shape=signTable$P.Val.Method, colour = logb(signTable$harmonicmean, 2),size = (signTable$Permutations)))+ 
        scale_colour_gradientn( colours = spetralc, breaks = legendbreaks, labels = legendlables, guide = "colourbar", limits = c(0, 9.2)) + ylim(0,60)+# xlim(0,xlimits[i])+
        scale_size_continuous(limits = c(10,10001),breaks = c(10,2000,5000,10001),labels = c("10","2000","5000","10000"))+
        geom_errorbar(aes(ymin=-log(P.value + CI.Upper, base = 2), ymax=-log(P.value - CI.Lower, base = 2)))+
        scale_shape_manual( values = c(2,3,4,1), 
                        breaks = c( "p_ecdf",  "p_ecdf_with_pseudo", "p_ecdf_with_pseudo (p_gpd=0)", "p_gpd" ),
                        labels = c( "Pecdf", "Pecdf with pseudo", "Pecdf with_pseudo (pgpd=0)", "Pgpd" ))+
        ggtitle(clades[i]) + ylab(TeX("$Surprisal (-log_2(p))$")) + xlab(xlabels[i]) + labs(color = "Sample Size", shape = "PmethodType",size = "Permutation number") + 
        theme(
          axis.text = element_text(size = 13, family = "serif"),
          axis.title = element_text(size = 14, family = "serif"),
          legend.text = element_text(size = 12.5, family = "serif"),
          legend.title = element_text(size = 14, family = "serif"),
          plot.title = element_text(hjust = 0.5, family = "serif", size = 15, face = "bold" ), 
          legend.key.height = unit(0.9, "cm"),
          plot.margin=unit(c(0.2,0.2,0.5,0.2), "cm")
        )
      plotls[[i]] <- p  
    }
    setEPS()
    postscript("KLD-Significance_ENRIETTII_MAJOR.eps", width = 16, height = 8,fonts=c("serif", "Palatino"))
    p <- ggarrange( plotls[[1]], plotls[[2]], plotls[[3]], plotls[[4]],ncol = 2, nrow = 2, common.legend = TRUE,legend.grob =         get_legend(plotls[[2]]), legend = "right",labels = c("A", "B", "C","D"))
    print(p)
    dev.off()
    ```



2. Creating the Supplemental Figure <ID-Significance_ENRIETTII_MAJOR.eps>

    1. Create the table of statistics for ID logos of clade `MAJOR` vs `HOMO` using the pvalue-calculation-method `GPD` with maximum permutation number = 10000, target permutation number 500, and calculate the Confidence Interval of p-values with 5% significance level. This will use the default entropy estimator NSB with default Maximum conditional sample size 5. The running time on an Intel Core i7 Dell XPS was ~4.76 hrs. The outputs are two text files: ID_MAJOR_HOMO_stats.txt and ID_HOMO_MAJOR_stats.txt
    ```shell
    tsfm --idperms 10000 -m GPD -c tRNA_L_skel_Leish.sites74.struct.cove HOMO/HOMO MAJOR/MAJOR
    mv ID_MAJOR_HOMO_stats.txt ID_HOMO_MAJOR_stats.txt GPD
    ```
    
    2. Create the table of statistics for ID logos of clade `ENRIETTII` vs `HOMO` using the pvalue-calculation-method `GPD` with maximum permutation number = 10000, target permutation number 500, and calculate the Confidence Interval of p-values with 5% significance level. Default -e NSB -x 5. The running time on an Intel Core i7 Dell XPS was ~4.89 hrs. The outputs are two text files: ID_ENRIETTII_HOMO_stats.txt and ID_HOMO_ENRIETTII_stats.txt
    ```shell
    tsfm --idperms 10000 -m GPD -c tRNA_L_skel_Leish.sites74.struct.cove HOMO/HOMO ENRIETTII/ENRIETTII
    mv ID_ENRIETTII_HOMO_stats.txt ID_HOMO_ENRIETTII_stats.txt GPD
    ```
    
    3. Create the figure using the ID tables in GPD folder by running this R script:
    ```R
    library(latex2exp)
    library(ggplot2)
    library(ggpubr)
    #setwd("PATH/TO/WHERE/THE/FOLDERS/GPD/ECDF/ECDF_pseudo/ARE/LOCATED") # make sure to set your working directory to where the input data folders GPD, ECDF and ECDF_pseudo are located.
    spetralc <-c("#D53E4F" ,"#F46D43" ,"#FDAE61" , "#66C2A5" ,"#3288BD","#313695")
    breaksarray <- c(1,2,4,10,30,80,200,500)
    legendbreaks <- logb(breaksarray, base = 2)
    legendlables <- as.character(breaksarray)
    ID_paths <-
      c(
        "GPD/ID_MAJOR_HOMO_stats.txt",
        "GPD/ID_HOMO_MAJOR_stats.txt",
        "GPD/ID_ENRIETTII_HOMO_stats.txt",
        "GPD/ID_HOMO_ENRIETTII_stats.txt"
      )
    xlabels <- c("Information Difference","Information Difference","Information Difference","Information Difference")
    clades <- c("MAJOR (n=664 sequences) against Homo (m=431 sequences)", "HOMO against MAJOR", "ENRIETTII (n=160 sequences) against Homo (m=431 sequences)", "HOMO against ENRIETTII")
    plotls <- list()
    for (i in 1:length(ID_paths)) {
      signTable <- read.table(ID_paths[i], header = TRUE, sep = "\t")
      signTable <- signTable[signTable$Sample.Sz.Back != 0 & signTable$Sample.Sz.Fore != 0 & signTable$Statistic != 0, ]
      signTable$harmonicmean <- 0
      for (j in 1:nrow(signTable)) {
        signTable$harmonicmean[j] = 1 / mean(1 / c(signTable$Sample.Sz.Back[j], signTable$Sample.Sz.Fore[j]))
      }
      p <- ggplot(signTable, aes_string( x = round(signTable$Statistic, 4), y = -log(signTable$P.value, base = 2),colour =         logb(signTable$harmonicmean, 2))) + 
        geom_point(aes_string( shape=signTable$P.Val.Method, colour = logb(signTable$harmonicmean, 2),size = (signTable$Permutations)))+ 
        scale_colour_gradientn( colours = spetralc, breaks = legendbreaks, labels = legendlables, guide = "colourbar", limits = c(0, 9.2)) + ylim(0,60)+# xlim(0,xlimits[i])+
        scale_size_continuous(limits = c(10,10001),breaks = c(10,2000,5000,10001),labels = c("10","2000","5000","10000"))+
        geom_errorbar(aes(ymin=-log(P.value + CI.Upper, base = 2), ymax=-log(P.value - CI.Lower, base = 2)))+
        scale_shape_manual( values = c(2,3,4,1), 
                        breaks = c( "p_ecdf",  "p_ecdf_with_pseudo", "p_ecdf_with_pseudo (p_gpd=0)", "p_gpd" ),
                        labels = c( "Pecdf", "Pecdf with pseudo", "Pecdf with_pseudo (pgpd=0)", "Pgpd" ))+
        ggtitle(clades[i]) + ylab(TeX("$Surprisal (-log_2(p))$")) + xlab(xlabels[i]) + labs(color = "Sample Size", shape = "PmethodType",size = "Permutation number") + 
        theme(
          axis.text = element_text(size = 13, family = "serif"),
          axis.title = element_text(size = 14, family = "serif"),
          legend.text = element_text(size = 12.5, family = "serif"),
          legend.title = element_text(size = 14, family = "serif"),
          plot.title = element_text(hjust = 0.5, family = "serif", size = 15, face = "bold" ), 
          legend.key.height = unit(0.9, "cm"),
          plot.margin=unit(c(0.2,0.2,0.5,0.2), "cm")
        )
      plotls[[i]] <- p  
    }
    setEPS()
    postscript("ID-Significance_ENRIETTII_MAJOR.eps", width = 16, height = 8,fonts=c("serif", "Palatino"))
    p <- ggarrange( plotls[[1]], plotls[[2]], plotls[[3]], plotls[[4]],ncol = 2, nrow = 2, common.legend = TRUE,legend.grob =         get_legend(plotls[[2]]), legend = "right",labels = c("A", "B", "C","D"))
    print(p)
    dev.off()
    ```

2. Creating the slope graph from tSFM manuscript. 

    1. Create the table of statistics for KLD logos of clade `ENRIETTII` vs `HOMO` using the pvalue-calculation-method `ECDF_pseudo` (naive Monte Carlo method using pseudo-counts) with maximum permutation number = 10000. The running time on a Linux compute node with 20 cores at 2301 MHz and 120 GB RAM was ~7.72 hrs. The outputs are two text files: KLD_ENRIETTII_HOMO_stats.txt and KLD_HOMO_ENRIETTII_stats.txt
    ```shell
    mkdir ECDF_pseudo
    tsfm --kldperms 10000 -m ECDF_pseudo -c tRNA_L_skel_Leish.sites74.struct.cove --clade HOMO HOMO/HOMO MAJOR/MAJOR
    mv KLD_ENRIETTII_HOMO_stats.txt KLD_HOMO_ENRIETTII_stats.txt ECDF_pseudo
    ```
    
    2. Create the table of statistics for KLD logos of clade `ENRIETTII` vs `HOMO` using the pvalue-calculation-method `ECDF` (Monte Carlo method terminating after M = 10 exceedances) with maximum permutation number = 10000. The running time on a Linux compute node with 20 cores at 2301 MHz and 120 GB RAM was ~4.86 hrs. The outputs are two text files: KLD_ENRIETTII_HOMO_stats.txt and KLD_HOMO_ENRIETTII_stats.txt
    ```shell
    mkdir ECDF
    tsfm --kldperms 10000 -m ECDF -c tRNA_L_skel_Leish.sites74.struct.cove --clade HOMO HOMO/HOMO MAJOR/MAJOR
    mv KLD_ENRIETTII_HOMO_stats.txt KLD_HOMO_ENRIETTII_stats.txt ECDF
    ```
    
    3. Create the figure using the KLD tables for ENRIETTII against HOMO from three folders GPD, ECDF and ECDF_pseudo
    ```R
    library(latex2exp)
    library(ggplot2)
    library(ggpubr)
    #setwd("PATH/TO/WHERE/THE/FOLDERS/GPD/ECDF/ECDF_pseudo/ARE/LOCATED") # make sure to set your working directory to where the input data folders GPD, ECDF and ECDF_pseudo are located.
    signTable0 <- read.table("ECDF_pseudo/KLD_ENRIETTII_HOMO_stats.txt", header = TRUE, sep = "\t")
    signTable1 <- read.table("ECDF/KLD_ENRIETTII_HOMO_stats.txt", header = TRUE, sep = "\t")
    signTable2 <- read.table("GPD/KLD_ENRIETTII_HOMO_stats.txt", header = TRUE, sep = "\t")
    signTable0 <- signTable0[signTable0$Sample.Sz.Back != 0 & signTable0$Sample.Sz.Fore != 0 & signTable0$Statistic != 0,]
    signTable0$method <- "ECDF-with-pseudo"
    signTable1 <- signTable1[signTable1$Sample.Sz.Back != 0 & signTable1$Sample.Sz.Fore != 0 & signTable1$Statistic != 0,]
    signTable1$method <- "ECDF"
    signTable2 <- signTable2[signTable2$Sample.Sz.Back != 0 & signTable2$Sample.Sz.Fore != 0 & signTable2$Statistic != 0,]
    signTable2$method <- "GPD with Initial-Target-Permnum = 500"
    signTable0$P.Val.Method <- "p_ecdf_with_pseudo" 
    common_cols <- intersect(colnames(signTable0), colnames(signTable1))
    bindeddf <- rbind(
      subset(signTable0, select = common_cols), 
      subset(signTable1, select = common_cols),
      subset(signTable2, select = common_cols)
    )
    bindeddf$coordstate <- paste(bindeddf$Coord,bindeddf$State,sep = "")
    bindeddf <- bindeddf[order(bindeddf$coordstate),]
    bindeddf$method <- factor(bindeddf$method,levels = c("ECDF-with-pseudo","ECDF","GPD with Initial-Target-Permnum = 500"))
    bindeddf$harmonicmean <- 0
    for (j in 1:nrow(bindeddf)) {
      bindeddf$harmonicmean[j] = 1 / mean(1 / c(bindeddf$Sample.Sz.Back[j], bindeddf$Sample.Sz.Fore[j]))
    }
    spetralc <-c("#D53E4F" ,"#F46D43" ,"#FDAE61" , "#66C2A5" ,"#3288BD","#313695")
    breaksarray <- c(1,2,4,10,30,80,200,500)
    legendbreaks <- logb(breaksarray, base = 2)
    legendlables <- as.character(breaksarray)
    p <- ggplot(data = bindeddf, aes( x = method ,  y =  jitter(-log(P.value, base = 2),factor = 80), group = coordstate )) +
      geom_line(aes(color = logb(harmonicmean, base = 2)), size = 0.5) +
      geom_point(aes(color = logb(harmonicmean, base = 2),shape=factor(P.Val.Method)), size = 2) +
      scale_x_discrete(position = "top") + 
      scale_colour_gradientn(  colours = spetralc, breaks = legendbreaks, labels = legendlables, guide = "colourbar", limits = c(0, 9.2) ) + 
      theme( legend.key.height = unit(0.7, 'in'), plot.title = element_text(hjust = 0.5), text = element_text(size = 24,  family = "serif")) + 
      scale_shape_manual(  values = c(2,3,4,1),  breaks = c( "p_ecdf", "p_ecdf_with_pseudo",  "p_ecdf_with_pseudo (p_gpd=0)", "p_gpd"),
                           labels = c( "pecdf", "pecdf with pseudo", "pecdf with pseudo (pgpd=0)", "pgpd" ))+
      labs(title = "" , colour = "Sample Size",shape="Pmethod", size = "Permutation number") +
      xlab("") + ylab(TeX("$Surprisal (-log_2(p))$"))
    setEPS()
    postscript("KLD_ENRIETTII_Slopegraph_GPDvsECDF.eps", width = 18, height = 10,fonts=c("serif", "Palatino"))
    print(p)
    dev.off()
    ```
 
