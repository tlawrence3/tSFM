# -*- coding: utf-8 -*-
import argparse
import sys
import os
import itertools
import tsfm.MolecularInformation as MolecularInformation
from tsfm._version import __version__

def main():
    # Setup parser
    parser = argparse.ArgumentParser(
        description="tSFM (tRNA Structure-Function Mapper) calculates functional Class-Informative Features (CIFs) and their evolutionary divergences for tRNAs or other RNA families.",
        epilog="Please cite Lawrence et al. (2020) tSFM: tRNA Structure-Function Mapper.")
    # Required arguments
    parser.add_argument("file_prefix",
                        help="One or more paths/file-prefix strings corresponding to sets of input files compiled for a single clade in clustalW format. Input files should be named <path>/<prefix>_<functional-class>.<extension>, where <functional_class> is a single letter",
                        nargs='+')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--infernal", type=argparse.FileType("r"),
                       help="Use secondary structure file INFERNAL, required to calculate functional information of base-pair features, in Infernal format")
    group.add_argument("-c", "--cove", type=argparse.FileType("r"),
                       help="Use secondary structure file COVE, required to calculate functional information of base-pair features, in COVE format. Example: \"#=CS  >>>>>>>..>>>>...........<<<<.>>>>>.......<<<<<.....>>>>>....\n#=CS      ...<<<<<<<<<<<<.\"")
    group.add_argument("-t", "--text", type=argparse.FileType("r"),
                       help="Use secondary structure file TEXT, required to calculate functional information of base-pair features, is in text format. Example: \"A:0,72,1,71,2,70,3,69,4,68,5,67,6,66\nD:9,25,10,24,11,23,12,22\nC:27,43,28,42,29,41,30,40,31,39\nT:49,65,50,64,51,63,52,62,53,61\"")
    group.add_argument("-s", "--single", action="store_true",
                       help="Do not calculate functional information for paired features. Calculate for single-site features only.")

    # group.add_argument("-f", "--file",     action="store_true", help="Read in previous results from file ")

    # Options
    parser.add_argument("-n", "--nosingle", action="store_true",
                       help="Do not calculate functional information for single-site features. Calculate for paired-site features only.")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s v{}".format(__version__))
    parser.add_argument("-p", "--processes", type=int, default=os.cpu_count(),
                        help="Set the maximum number of concurrent processes. Default is the number of cores reported by the operating system.")
    parser.add_argument("-e", "--entropy", type=str, default="NSB",
                        help="Use entropy estimator ENTROPY when conditional sample sizes exceed maximum for exact calculation. If value is \"NSB\", use Nemenman-Shafee-Bialek estimator. If value is \"MM\", use Miller-Madow estimator. Default is NSB",
                        choices=['NSB', 'MM'])
    parser.add_argument("-x", "--exact",
                        help="Maximum conditional sample size to exactly calculate entropy correction. Default = 5",
                        type=int, default=5)
    parser.add_argument("-l", "--logos",
                        help="Visualize feature information using function/inverse function logos in extended postscript format.",
                        action="store_true")
    parser.add_argument("-v", "--inverse", action="store_true",
                        help="Additionally calculate functional information for inverse features/logos (features under-represented in specific functional classes)")
    parser.add_argument("-P", "--permutations",
                        help="Calculate the significance of CIFs by a permutation test, with a number of permutations equal to PERMUTATIONS (an integer). Default is to not calculate significance of CIFs.",
                        type=int, default=0)
    parser.add_argument("-C",
                        help="Specify a method for multiple test correction for significance calculations: bonferroni, sidak, holm, holm-sidak, simes-hochberg, hommel, BH (Benjamini-Hochberg FDR), BY (Benjamini-Yekutieli FDR) or GBS (Gavrilov-Benjamini-Sarkar FDR). Default is BH", default="BH", choices=['bonferroni', 'sidak', 'holm', 'holm-sidak', 'simes-hochberg', 'hommel', 'BH', 'BY','GBS'], dest="correction")   
    parser.add_argument("-T",
                        help="Test the significance of only CIF stack-heights, only CIF letter-heights, or both. Default is both.",
                        default="both",
                        choices=['stacks', 'letters', 'both'], dest="test") 
    parser.add_argument("-I", "--idlogos", help='Compute Information Differece logos for each pair of clades',
                        action="store_true")
    parser.add_argument("-K", "--kldlogos", help='Compute Kullback-Liebler Divergence logos for each pair of clades',
                        action="store_true")
    parser.add_argument("-B", "--bubbles",
                        help='Compute input table for structural bubble-plots (to be computed in R) like those appearing in Kelly et al. (2020).',
                        action="store_true")
    parser.add_argument("--kldperms",
                        help="Set the number of permutations to compute significance of Kullback-Leibler Divergences. Default is to not calculate signifiance.",
                        type=int, default=0)
    parser.add_argument("--idperms",
                        help="Set the number of permutations to compute significance of Information Differences (SLOW). Default is to not calculate signifiance.",
                        type=int, default=0)
    parser.add_argument("-J", "--JSD",
                        help="Produce pairwise distance matrices between function logos for different taxa based on Jensen-Shannon Divergences",
                        action="store_true")
    parser.add_argument("--clade", type=str,
                        help="Contrast clade CLADE against all others. CLADE should be one of the file-prefix strings passed as a required argument to the program, stripped of its path. Default is to compute contrasts for all pairs of clades.")

    args = parser.parse_args()

    if (args.single and args.nosingle):
        sys.exit("Options --single and --nosingle are incompatible.")
    
    # initialize dictionary that contains all datasets labeled by the file prefix
    logo_dict = {}

    # create results object from previously calculated function logos
    # if (args.file):
    #    results = {}
    #    for prefix in args.file_prefix:
    #        prefix_name = prefix.split("/")[-1]
    #        results[prefix_name] = MolecularInformation.FunctionLogoResults(prefix, from_file = True)

    # load datasets into logo objects
    # else:

    if (args.text):
        for prefix in args.file_prefix:
            prefix_name = prefix.split("/")[-1]
            logo_dict[prefix_name] = MolecularInformation.FunctionLogo(args.text, "text", None, None)
    elif (args.cove):
        for prefix in args.file_prefix:
            prefix_name = prefix.split("/")[-1]
            logo_dict[prefix_name] = MolecularInformation.FunctionLogo(args.cove, "cove", None, None)
    elif (args.infernal):
        for prefix in args.file_prefix:
            prefix_name = prefix.split("/")[-1]
            logo_dict[prefix_name] = MolecularInformation.FunctionLogo(args.infernal, "infernal", None, None)
    elif (args.single):
        for prefix in args.file_prefix:
            prefix_name = prefix.split("/")[-1]
            logo_dict[prefix_name] = MolecularInformation.FunctionLogo(None, "s", None, None)

    for prefix in args.file_prefix:
        prefix_name = prefix.split("/")[-1]
        logo_dict[prefix_name].parse_sequences(prefix)

    if (args.clade and args.clade not in logo_dict.keys()):
        sys.exit("tsfm: Argument to option --clade must be identical to one of the file-prefix arguments to the program, stripped of its path.")

    if (args.bubbles and not args.clade):
        sys.exit("tsfm: Option --bubbles requires designation of a specific clade to contrast against using option --clade.")
        
        

    # Calculate exact method sample size correction
    if (args.exact):
        for key in logo_dict:
            print("Calculating Sample Size Correction for {}".format(key))
            logo_dict[key].calculate_exact(args.exact, args.processes)
            if (args.inverse):
                print("Calculating Sample Size Correction for Inverse {}".format(key))
                logo_dict[key].calculate_exact(args.exact, args.processes, inverse=True)

    # Perform function label swapping permutations and calculate entropy distribution from permutations
    if (args.permutations):
        multitest_methods = {'bonferroni': 'b', 'sidak': 's', 'holm': 'h',
                             'holm-sidak': 'hs', 'simes-hochberg': 'sh',
                             'hommel': 'ho', 'BH': 'fdr_bh', 'BY': 'fdr_by',
                             'GBS': 'fdr_gbs'}
        perm_dict = {}
        for key in logo_dict:
            print("Generating permuted alignment data for {}".format(key), file=sys.stderr)
            logo_dict[key].permute(args.permutations, args.processes)
        for key in logo_dict:
            print("Calculating permutation information for {}".format(key), file=sys.stderr)
            perm_dict[key] = logo_dict[key].permInfo(args.entropy, args.processes)
        if (args.inverse):
            perm_inverse_dict = {}
            for key in logo_dict:
                print("Calculating inverse permutation information for {}".format(key), file = sys.stderr)
                perm_inverse_dict[key] = logo_dict[key].permInfo(args.entropy, args.processes, inverse=True)

    results = {}

    # Initialization of function logo result objects
    for key in logo_dict:
        results[key] = MolecularInformation.FunctionLogoResults(key,
                                                                logo_dict[key].basepairs,
                                                                logo_dict[key].pos,
                                                                logo_dict[key].sequences,
                                                                logo_dict[key].pairs,
                                                                logo_dict[key].singles)

    if (args.entropy == "NSB"):
        for key in logo_dict:
            print("Calculating information statistics for {} using NSB estimator".format(key), file = sys.stderr)
            info, height_dict = logo_dict[key].calculate_entropy_NSB()
            results[key].add_information(info=info, height=height_dict)
            if (args.inverse):
                print("Calculating inverse information statistics for {} using NSB estimator".format(key), file
                      = sys.stderr)
                info_inverse, height_dict_inverse = logo_dict[key].calculate_entropy_inverse_NSB()
                results[key].add_information(info=info_inverse, height=height_dict_inverse, inverse=True)
    elif (args.entropy == "MM"):
        for key in logo_dict:
            print("Calculating information statistics using Miller-Maddow estimator")
            info, height_dict = logo_dict[key].calculate_entropy_MM()
            results[key].add_information(info=info, height=height_dict)
            if (args.inverse):
                print("Calculating inverse using Miller-Maddow estimator")
                info_inverse, height_dict_inverse = logo_dict[key].calculate_entropy_inverse_MM()
                results[key].add_information(info=info_inverse, height=height_dict_inverse, inverse=True)

    if (args.permutations):
        print("Calculating p-values using {} multiple test correction".format(args.correction))
        for key in results:
            results[key].add_stats(perm_dict[key], multitest_methods[args.correction], args.test, args.nosingle)
            if (args.inverse):
                results[key].add_stats(perm_inverse_dict[key], multitest_methods[args.correction], args.test,
                                       args.nosingle, inverse=True)

    for key in results:
        print("Writing text output for {}".format(key))
        results[key].text_output(args.correction)

    if (args.logos and not args.inverse):
        for key in results:
            print("Writing function logo postscript files for {}".format(key))
            results[key].logo_output(logo_prefix="functionlogo")
    elif (args.logos and args.inverse):
        for key in results:
            print("Writing inverse function logo postscript files for {}".format(key))
            results[key].logo_output(inverse=True)

    if (args.JSD):
        distance = MolecularInformation.DistanceCalculator("jsd")
        distance.get_distance(results)

    # ______________________________________________________________________________________________________________
    if args.kldlogos or args.idlogos or args.bubbles or args.idperms or args.kldperms:
        info_height_dic = {}
        for key in results:
            info_height_dic[key] = {"info": results[key].info, "height": results[key].height}

        # Initializing variables of class FunctionLogoDifference  __________________________________________________
        pos = results[list(logo_dict.keys())[0]].pos
        basepair = results[list(logo_dict.keys())[0]].basepairs
        types = logo_dict[list(logo_dict.keys())[0]].functions
        # Variables pairs and single will be initialized later for each pair of clades separately __________________

        results_prob_dist = {}
        post_nopseudo = {}
        if args.clade:
            pairwise_combinations = [(x, y) for (x, y) in itertools.product(
                [element for element in logo_dict.keys() if element not in args.clade], [args.clade])]
        else:
            pairwise_combinations = itertools.combinations(logo_dict.keys(), 2)

        for cpair in pairwise_combinations:
            print("Calculating ID and/or KLD for", cpair[0], "and", cpair[1])
            pairwise_permutation = itertools.permutations(list(cpair), 2)
            for pair in pairwise_permutation:
                pairs = list(set(logo_dict[pair[0]].pairs) & set(logo_dict[pair[1]].pairs))
                single = list(set(logo_dict[pair[0]].singles) & set(logo_dict[pair[1]].singles))
                difference = MolecularInformation.FunctionLogoDifference(pos, types, pairs, basepair, single)

                results_prob_dist[pair[0]] = {}
                results_prob_dist[pair[0]]['post'], results_prob_dist[pair[0]][
                    'prior'] = difference.calculate_prob_dist_pseudocounts(logo_dict[pair[0]], logo_dict[pair[1]])
                post_nopseudo[pair[0]] = difference.calculate_prob_dist_nopseudocounts(logo_dict[pair[0]])
            kld_height_dic = {}  # KLDs are saved with the background key
            ratios_dic = {}  # ratios are saved with the background key
            id_height_dic = {}  # IDs are saved with background key
            kld_infos = {}
            id_infos = {}
            pairwise_permutation = itertools.permutations(list(cpair), 2)
            for pair in pairwise_permutation:
                pairs = list(set(logo_dict[pair[0]].pairs) & set(logo_dict[pair[1]].pairs))
                single = list(set(logo_dict[pair[0]].singles) & set(logo_dict[pair[1]].singles))
                difference = MolecularInformation.FunctionLogoDifference(pos, types, pairs, basepair, single)

                ratios_dic[pair[0]] = difference.calculate_ratios(back_prior=results_prob_dist[pair[0]]['prior'],
                                                                  fore_prior=results_prob_dist[pair[1]]['prior'],
                                                                  back_post=results_prob_dist[pair[0]]['post'],
                                                                  nopseudo_post_fore=post_nopseudo[pair[1]])
                if args.kldlogos or args.bubbles or args.kldperms:
                    kld_info, kld_height = difference.calculate_kld(logo_dict, key_back=pair[0], key_fore=pair[1],
                                                                    back_prior=results_prob_dist[pair[0]]['prior'],
                                                                    fore_prior=results_prob_dist[pair[1]]['prior'],
                                                                    back_post=results_prob_dist[pair[0]]['post'],
                                                                    fore_post=results_prob_dist[pair[1]]['post'],
                                                                    ratios=ratios_dic[pair[0]])
                    if args.kldlogos:
                        print("Writing KLD logos for", cpair[0], "and", cpair[1])
                        logoprefix = "KLDlogo"
                        results[pair[0]].add_information(info=kld_info, height=kld_height)
                        results[pair[0]].logo_output(logo_prefix=logoprefix, logo_postfix=pair[1])

                    kld_height_dic[pair[0]] = {"kld": kld_info, "height": kld_height}
                    kld_infos[pair[0]] = kld_info

                if args.idlogos or args.bubbles or args.idperms:
                    id_info = difference.calculate_logoID_infos(info_b=info_height_dic[pair[0]]['info'],
                                                                info_f=info_height_dic[pair[1]]['info'])
                    id_height = difference.calculate_logoID_heights(info=id_info, ratios=ratios_dic[pair[0]])
                    id_height_dic[pair[0]] = {"id": id_info, "height": id_height}

                    if args.idlogos:
                        print("Writing ID logos for", cpair[0], "and", cpair[1])
                        results[pair[0]].add_information(info=id_info, height=id_height)
                        logoprefix = "IDlogo"
                        results[pair[0]].logo_output(logo_prefix=logoprefix, logo_postfix=pair[1])
                    id_infos[pair[0]] = id_info

            if args.bubbles:
                pairwise_permutation = itertools.permutations(list(cpair), 2)
                print("Writing text output for bubble plots")
                for pair in pairwise_permutation:
                    # pair[0] is background
                    pairs = list(set(logo_dict[pair[0]].pairs) & set(logo_dict[pair[1]].pairs))
                    single = list(set(logo_dict[pair[0]].singles) & set(logo_dict[pair[1]].singles))
                    difference = MolecularInformation.FunctionLogoDifference(pos, types, pairs, basepair, single)

                    difference.func_ID_KLD_2table(fore_logo_info=info_height_dic[pair[1]]['info'],
                                                  fore_logo_height=info_height_dic[pair[1]]['height'],
                                                  fore_idlogo_info=id_height_dic[pair[0]]['id'],
                                                  back_idlogo_info=id_height_dic[pair[1]]['id'],
                                                  fore_idlogo_height=id_height_dic[pair[0]]['height'],
                                                  back_idlogo_height=id_height_dic[pair[1]]['height'],
                                                  kld_info=kld_height_dic[pair[0]]['kld'],
                                                  kld_height=kld_height_dic[pair[0]]['height'],
                                                  fore=pair[1],back=pair[0])
            if args.kldperms:
                print("Calculating significance of KLDs between", cpair[0], "and", cpair[1])
                pairs = list(set(logo_dict[cpair[0]].pairs) & set(logo_dict[cpair[1]].pairs))
                single = list(set(logo_dict[cpair[0]].singles) & set(logo_dict[cpair[1]].singles))
                klddifference = MolecularInformation.FunctionLogoDifference(pos, types, pairs, basepair, single)
                logo_dict_pair = {key: logo_dict[key] for key in [cpair[0], cpair[1]]}
                kld_pvalues = klddifference.calculate_kld_significance(logo_dict_pair, kld_infos, args.kldperms,
                                                                       args.processes)
                kld_pvalues_corrected = klddifference.addstats(kld_pvalues, args.correction)

                print("Writing text output for KLD significance")
                klddifference.write_pvalues(kld_pvalues, kld_pvalues_corrected, kld_infos, logo_dict_pair, "KLD")

            if args.idperms:
                    print("Calculating significance of IDs between", cpair[0], "and", cpair[1])
                    pairs = list(set(logo_dict[cpair[0]].pairs) & set(logo_dict[cpair[1]].pairs))
                    single = list(set(logo_dict[cpair[0]].singles) & set(logo_dict[cpair[1]].singles))
                    iddifference = MolecularInformation.FunctionLogoDifference(pos, types, pairs, basepair, single)
                    logo_dict_pair = {key: logo_dict[key] for key in [cpair[0], cpair[1]]}
                    id_pvalues = iddifference.calculate_id_significance(logo_dict_pair, id_infos, args.idperms,
                                                                        args.processes,
                                                                        args.exact,
                                                                        args.entropy)
                    id_pvalues_corrected = iddifference.addstats(id_pvalues, args.correction)
                    print("Writing text output for ID significance")
                    iddifference.write_pvalues(id_pvalues, id_pvalues_corrected, id_infos, logo_dict_pair, "ID")

if __name__ == "__main__":
    main()
