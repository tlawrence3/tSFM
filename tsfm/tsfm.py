# -*- coding: utf-8 -*-
import argparse
import sys
import os
import tsfm.MolecularInformation as MolecularInformation
from tsfm._version import __version__

def main():
    # Setup parser
    parser = argparse.ArgumentParser(description="tsfm")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s v{}".format(__version__))
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--infernal", type=argparse.FileType("r"), help="Structure file is in infernal format")
    group.add_argument("-c", "--cove", type=argparse.FileType("r"), help="Structure file is in cove format")
    group.add_argument("-t", "--text", type=argparse.FileType("r"), help="Structure file is in text format")
    group.add_argument("-f", "--file", action="store_true", help="use to read in previous results from file")
    group.add_argument("-s", "--single", action="store_true",
                       help="Calculate functional information for single sites only")
    parser.add_argument("-p", "--proc", type=int, default=os.cpu_count(),
                        help="Maximum number of concurrent processes. Default is number of cores reported by the OS")
    parser.add_argument("-v", "--inverse", action="store_true", help="calculate anti-determinates")
    parser.add_argument("-a", "--alpha", type=float, default=0.05, help="Current not implemented. Default = 0.05")
    parser.add_argument("-e", "--entropy", type=str, default="NSB",
                        help="Method of entropy estimation. Either NSB or Miller. Default = NSB")
    parser.add_argument('--max', '-x',
                        help="Maximum sample size to calculate the exact entropy correction. Default = 10", type=int,
                        default=10)
    parser.add_argument('--logo', help='Produce function logo ps files', action="store_true")
    parser.add_argument("-B", help="Number of permutations. Default value is 100", type=int, default=0)
    parser.add_argument("-M",
                        help="Specify method to correct p-values for multiple-comparisons. Current methods available: bonferroni, sidak, holm, holm-sidak, simes-hochberg, hommel, BH, BY, TSBH, TSBKY, and GBS. Default is BH",
                        default="BH")
    parser.add_argument("-j", "--jsd", action="store_true", help="")
    parser.add_argument("file_prefix", help="File prefix", nargs='+')
    parser.add_argument("--idlogo", help='id logo', action="store_true")
    parser.add_argument("--kldlogo", help='KLD logo', action="store_true")
    parser.add_argument("--bt", help='bubble table', action="store_true")
    args = parser.parse_args()

    logo_dict = {}

    if (args.file):
        results = {}
        for prefix in args.file_prefix:
            prefix_name = prefix.split("/")[-1]
            results[prefix_name] = MolecularInformation.FunctionLogoResults(prefix, from_file=True)
    else:
        if (args.text):
            for prefix in args.file_prefix:
                prefix_name = prefix.split("/")[-1]
                logo_dict[prefix_name] = MolecularInformation.FunctionLogo(args.text, "text")
        elif (args.cove):
            for prefix in args.file_prefix:
                prefix_name = prefix.split("/")[-1]
                logo_dict[prefix_name] = MolecularInformation.FunctionLogo(args.cove, "cove")
        elif (args.single):
            for prefix in args.file_prefix:
                prefix_name = prefix.split("/")[-1]
                logo_dict[prefix_name] = MolecularInformation.FunctionLogo(args.cove, "s")

        for prefix in args.file_prefix:
            prefix_name = prefix.split("/")[-1]
            logo_dict[prefix_name].parse_sequences(prefix)

        if (args.max):
            for key in logo_dict:
                print("Calculating Sample Size Correction for {}".format(key))
                logo_dict[key].calculate_exact(args.max, args.proc)
                if (args.inverse):
                    print("Calculating Sample Size Correction for Inverse {}".format(key))
                    logo_dict[key].calculate_exact(args.max, args.proc, inverse=True)

        if (args.B):
            multitest_methods = {'bonferroni': 'b', 'sidak': 's', 'holm': 'h',
                                 'holm-sidak': 'hs', 'simes-hochberg': 'sh',
                                 'hommel': 'ho', 'BH': 'fdr_bh', 'BY': 'fdr_by',
                                 'TSHB': 'fdr_tsbh', 'TSBKY': 'fdr_tsbky',
                                 'GBS': 'fdr_gbs'}
            perm_dict = {}
            for key in logo_dict:
                print("Generating permuted alignment data for {}".format(key), file=sys.stderr)
                logo_dict[key].permute(args.B, args.proc)
            for key in logo_dict:
                print("Calculating permutation information for {}".format(key), file=sys.stderr)
                perm_dict[key] = logo_dict[key].permInfo(args.entropy, args.proc)
            if (args.inverse):
                perm_inverse_dict = {}
                for key in logo_dict:
                    print("Calculating inverse permutation information for {}".format(key), file=sys.stderr)
                    perm_inverse_dict[key] = logo_dict[key].permInfo(args.entropy, args.proc, inverse=True)

        results = {}
        for key in logo_dict:
            results[key] = MolecularInformation.FunctionLogoResults(key,
                                                                    logo_dict[key].basepairs,
                                                                    logo_dict[key].pos,
                                                                    logo_dict[key].sequences,
                                                                    logo_dict[key].pairs,
                                                                    logo_dict[key].singles)

        if (args.entropy == "NSB"):
            for key in logo_dict:
                print("Calculating information statistics for {} using NSB estimator".format(key), file=sys.stderr)
                info, height_dict = logo_dict[key].calculate_entropy_NSB()
                results[key].add_information(info=info, height=height_dict)
                if (args.inverse):
                    print("Calculating inverse information statistics for {} using NSB estimator".format(key),
                          file=sys.stderr)
                    info_inverse, height_dict_inverse = logo_dict[key].calculate_entropy_inverse_NSB()
                    results[key].add_information(info=info_inverse, height=height_dict_inverse, inverse=True)
        else:
            for key in logo_dict:
                print("Calculating information statistics using Miller-Maddow estimator")
                info, height_dict = logo_dict[key].calculate_entropy_MM()
                results[key].add_information(info=info, height=height_dict)
                if (args.inverse):
                    print("Calculating inverse using Miller-Maddow estimator")
                    info_inverse, height_dict_inverse = logo_dict[key].calculate_entropy_inverse_MM()
                    results[key].add_information(info=info_inverse, height=height_dict_inverse, inverse=True)

        if (args.B):
            print("Calculating p-values")
            for key in results:
                results[key].add_stats(perm_dict[key], multitest_methods[args.M])
                if (args.inverse):
                    results[key].add_stats(perm_inverse_dict[key], multitest_methods[args.M], inverse=True)

        for key in results:
            print("Writing text output for {}".format(key))
            results[key].text_output()

        if (args.logo and not args.inverse):
            for key in results:
                print("Writing function logo postscript files for {}".format(key))
                results[key].logo_output()
        elif (args.logo and args.inverse):
            for key in results:
                print("Writing function logo postscript files for {}".format(key))
                results[key].logo_output(inverse=True)

    if (args.jsd):
        distance = MolecularInformation.DistanceCalculator("jsd")
        distance.get_distance(results)
    # ____________________________________________________________________________________________________________
    # new id

    if args.idlogo or args.bt:
        distance = MolecularInformation.DistanceCalculator("id")
        heights_id, infos_id = distance.get_distance(results)  # heights are saved with their foreground key
    # ____________________________________________________________________________________________________________
    # new kld

    if args.kldlogo or args.bt:

        # ________________________

        # info_height_dic: info and heights for two keys
        info_height_dic = {}
        for key in results:
            info_height_dic[key] = {"info": results[key].info, "height": results[key].height}

        key_1 = list(logo_dict.keys())[0]
        key_2 = list(logo_dict.keys())[1]

        # validating variables:
        # single: which is based on the alphabets used in the trna sequences of each key.
        # Using the intersection in case one of the sequences has an unknown letter such as N.
        # pos: length of sequences which is similar for both keys (~= 72)
        # pairs: all the possible basepairs of alphabet. Using the intersection for the same reason as variable single
        # types: 21 functiona classes, assuming both organisms have at least one sequence for each class

        single = list(set(logo_dict[key_2].singles) & set(logo_dict[key_1].singles))
        pos = results[key_1].pos
        pairs = list(set(logo_dict[key_2].pairs) & set(logo_dict[key_1].pairs))
        basepair = results[key_1].basepairs
        types = logo_dict[key].functions

        # ________________________

        resultsKLD = {}  # dictionary of posterior and prior probability for each key
        for key in logo_dict.keys():
            prior, post = logo_dict[key].calculate_prob_dist(types, pairs=pairs, basepairs=basepair, singles=single)
            resultsKLD[key] = {}
            resultsKLD[key]['post'] = post
            resultsKLD[key]['prior'] = prior

        kld_height_dic = {}  # kld infos are saved with their foreground key same as ID infos
        pairwise_combinations = itertools.permutations(logo_dict.keys(), 2)
        for pair in pairwise_combinations:
            # pair[1] as foreground and pair[0] as background
            # exact is calculated for foreground in KLD So, we need to call calculate_kld with foreground object
            # so that we have that information for calculating exact

            kld_info, kld_height = logo_dict[pair[1]].calculate_kld(types, back_prior=resultsKLD[pair[0]]['prior'],
                                                                    fore_prior=resultsKLD[pair[1]]['prior'],
                                                                    back_post=resultsKLD[pair[0]]['post'],
                                                                    fore_post=resultsKLD[pair[1]]['post'], pairs=pairs,
                                                                    basepairs=basepair, singles=single)
            kld_height_dic[pair[1]] = {"kld": kld_info, "height": kld_height}

            # don't remember why callsed these two function with background key! does it matter?!
            results[pair[0]].add_information(info=kld_info, height=kld_height)
            results[pair[0]].logo_output()
        # bubble tables
        if args.bt:
            pairwise_combinations = itertools.permutations(logo_dict.keys(), 2)
            for pair in pairwise_combinations:
                # pair[1] as foreground and pair[0] as background
                logo_dict[pair[1]].func_ID_KLD_2table(types, fore_logo_info=info_height_dic[pair[1]]['info'],
                                                      fore_logo_height=info_height_dic[pair[1]]['height'],
                                                      fore_idlogo_info=infos_id[pair[1]],
                                                      back_idlogo_info=infos_id[pair[0]],
                                                      fore_idlogo_height=heights_id[pair[1]],
                                                      back_idlogo_height=heights_id[pair[0]],
                                                      kld_info=kld_height_dic[pair[1]]['kld'],
                                                      kld_height=kld_height_dic[pair[1]]['height'], functions=types,
                                                      states=single, pos=pos,
                                                      fore=pair[1])


if __name__ == "__main__":
    main()
