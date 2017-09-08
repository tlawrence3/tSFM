Tutorial
========

.. program:: tsfm

.. option:: -i (infernal), --infernal (infernal)

    Secondary structure annotation file in infernal format

.. option:: -c (cove), --cove (cove)

    Secondary structure annotation file in cove format

.. option:: -t (text), --text (text)

    Secondary structure annotation file in text format

.. option:: -f, --file

    Read previous results from file(s)

.. option:: -p <proc>. --proc <proc>

    Number of processor cores to use for concurrent processes. Default is
    the number of cores reported by the OS.

.. option:: -v, --inverse

    Calculate anti-determinates

.. option:: -a, --alpha

    Alpha value used to determine statistical significance for
    visualization and data filtering (not implemented).

.. option:: -e <entropy>, --entropy <entropy>

    Method of entropy estimation. NSB <NSB> and Miller-Maddow <Miller> 
    entropy esitmators are implemented. Default is NSB.

.. option:: -x <N>, --max <N>

    Calculate the exact method of small sample size correction for up to N samples. 
    This method is fully described in Schneider et al 1986. This calculation is
    polynomial in sample size. It becomes prohibitively expensive to calculate beyond
    a sample size of 16. The correction factor of each sample size will be calculated
    in parallel up to :option:`--proc` at a time. Default: 10.
 
.. option:: --logo

    Produce function logo postscript files

.. option:: -B <B>

    Number of permutations used to estimate discrete probability
    distributions for statisical testing. Default: 100.

.. option:: -M <M>

    Method to correct p-values for multiple-comparisons. Current 
    methods available: bonferroni, holm, hommel, BH, BY, and hochberg-simes. 
    Default value: BH.

.. option:: -j, --jsd

    Calculate the pairwise square root of the Jensen-Shannon divergence for all
    pairwise combinations of datasets.

.. option:: file_prefix

   Positional argument listing the file prefixes of the sequence alignment 
   files. Sequence alignment files are required to be in clustal format with 
   each functional class having its own file for each prefix. Alignment files 
   must conform to the naming standard fileprefix_functionalclass.aln.

   Example::

       To calculate functional information for the two datasets below 
       file_prefix would be: DataSet1 DataSet2

       - DataSet1_A.aln
       - DataSet1_C.aln
       - DataSet1_D.aln
       - DataSet1_E.aln
       - DataSet1_F.aln
       - DataSet1_G.aln
       - DataSet2_A.aln
       - DataSet2_C.aln
       - DataSet2_D.aln
       - DataSet2_E.aln
       - DataSet2_F.aln
       - DataSet2_G.aln
