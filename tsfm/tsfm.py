# -*- coding: utf-8 -*-
import argparse
import sys
import os
import itertools
import tsfm.MolecularInformation as MolecularInformation
from tsfm._version import __version__

def main():
     #Setup parser
    parser = argparse.ArgumentParser(description = "tsfm")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s v{}".format(__version__))
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--infernal", type=argparse.FileType("r"), help="Structure file is in infernal format")
    group.add_argument("-c", "--cove", type=argparse.FileType("r"), help="Structure file is in cove format")
    group.add_argument("-t", "--text", type=argparse.FileType("r"), help="Structure file is in text format")
    group.add_argument("-f", "--file", action="store_true", help="use to read in previous results from file")
    group.add_argument("-s", "--single", action="store_true", help="Calculate functional information for single sites only")
    parser.add_argument("-p", "--proc", type = int, default = os.cpu_count(), help="Maximum number of concurrent processes. Default is number of cores reported by the OS")
    parser.add_argument("-v","--inverse", action="store_true", help="calculate anti-determinates")
    parser.add_argument("-a", "--alpha", type=float, default=0.05, help="Current not implemented. Default = 0.05")
    parser.add_argument("-e", "--entropy", type=str, default = "NSB", help= "Method of entropy estimation. Either NSB or Miller. Default = NSB")
    parser.add_argument('--max', '-x', help="Maximum sample size to calculate the exact entropy correction. Default = 10", type=int, default = 0)
    parser.add_argument('--logo', help='Produce function logo ps files', action="store_true")
    parser.add_argument("-P", help="Number of permutations. Default value is 0", type=int, default=0)
    #Added command line argument for bootstrap reps
    parser.add_argument("-b", "--bootstrap", type=int, help="Number of bootstrap reps to perform. Default = 0", default = 0)
    parser.add_argument("-M", help = "Specify method to correct p-values for multiple-comparisons. Current methods available: bonferroni, sidak, holm, holm-sidak, simes-hochberg, hommel, BH, BY, TSBH, TSBKY, and GBS. Default is BH", default = "BH")
    parser.add_argument("-j", "--jsd", action="store_true", help="Produce jsd pairwise distance matrix between function logos")
    parse.add_argument("-I", "--ID", action="store_true", help="")
    parser.add_argument("-K", "--KLD", action="store_true", help="")
    parser.add_argument("file_prefix", help="File prefix", nargs='+')
    args = parser.parse_args()

    #dictionary that contains all datasets labeled by the file prefix
    logo_dict = {}

    #create results object from previously calculated function logos
    if (args.file):
        results = {}
        for prefix in args.file_prefix:
            prefix_name = prefix.split("/")[-1]
            results[prefix_name] = MolecularInformation.FunctionLogoResults(prefix, from_file = True)
    #load datasets into logo objects        
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

        #calculate exact method sample size correction
        if (args.max):
            for key in logo_dict:
                print("Calculating Sample Size Correction for {}".format(key))
                logo_dict[key].calculate_exact(args.max, args.proc)
                if (args.inverse):
                    print("Calculating Sample Size Correction for Inverse {}".format(key))
                    logo_dict[key].calculate_exact(args.max, args.proc, inverse=True)

        #Entry point for performing bootstrap replicates
        if (args.bootstrap):
            print("Generating bootstrap samples")
            for key in logo_dict:
                logo_dict[key].bootstrap(args.bootstrap, args.proc)

        #Perform function label swapping permutations and calculate entropy distribution from permuations
        if (args.P):
            multitest_methods = {'bonferroni': 'b', 'sidak': 's', 'holm': 'h',
                                 'holm-sidak': 'hs', 'simes-hochberg': 'sh',
                                 'hommel': 'ho', 'BH': 'fdr_bh', 'BY': 'fdr_by',
                                 'TSHB': 'fdr_tsbh', 'TSBKY': 'fdr_tsbky',
                                 'GBS': 'fdr_gbs'}
            perm_dict = {}
            for key in logo_dict:
                print("Generating permuted alignment data for {}".format(key), file=sys.stderr)
                logo_dict[key].permute(args.P, args.proc)
            for key in logo_dict:
                print("Calculating permutation information for {}".format(key), file=sys.stderr)
                perm_dict[key] = logo_dict[key].permInfo(args.entropy, args.proc)
            if (args.inverse):
                perm_inverse_dict = {}
                for key in logo_dict:
                    print("Calculating inverse permutation information for {}".format(key), file = sys.stderr)
                    perm_inverse_dict[key] = logo_dict[key].permInfo(args.entropy, args.proc, inverse = True)

        results = {}
        #Initialization function logo result objects
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
                results[key].add_information(info = info, height = height_dict)
                if (args.inverse):
                    print("Calculating inverse information statistics for {} using NSB estimator".format(key), file = sys.stderr)
                    info_inverse, height_dict_inverse = logo_dict[key].calculate_entropy_inverse_NSB()
                    results[key].add_information(info = info_inverse, height = height_dict_inverse, inverse = True)
        else:
            for key in logo_dict:
                print("Calculating information statistics using Miller-Maddow estimator")
                info, height_dict = logo_dict[key].calculate_entropy_MM()
                results[key].add_information(info = info, height = height_dict)
                if (args.inverse):
                    print("Calculating inverse using Miller-Maddow estimator")
                    info_inverse, height_dict_inverse = logo_dict[key].calculate_entropy_inverse_MM()
                    results[key].add_information(info = info_inverse, height = height_dict_inverse, inverse = True)

        if (args.P):
            print("Calculating p-values")
            for key in results:
                results[key].add_stats(perm_dict[key], multitest_methods[args.M])
                if (args.inverse):
                    results[key].add_stats(perm_inverse_dict[key], multitest_methods[args.M], inverse = True)

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
    if (args.ID):
        pairwise_combinations = itertools.permutations(results.keys(), 2)
        for pair in pairwise_combinations:
            pass
if __name__ == "__main__":
    main()
