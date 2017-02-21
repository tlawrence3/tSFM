# -*- coding: utf-8 -*-
import argparse
import os
import bplogofuntest.seq as seq

def main():
     #Setup parser
    parser = argparse.ArgumentParser(description = "bpLogoFun")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--infernal", action="store_true", help="Structure file is in infernal format")
    group.add_argument("-c", "--cove", action="store_true", help="Structure file is in cove format")
    group.add_argument("-t", "--text", action="store_true", help="Structure file is in text format")
    parser.add_argument("-p", "--proc", type = int, default = os.cpu_count(), help="Maximum number of concurrent processes. Default is number of cores reported by the OS")
    parser.add_argument("-v","--inverse", action="store_true", help="calculate anti-determinates")
    parser.add_argument("-a", "--alpha", type=float, default=0.05, help="Alpha value used for statistical tests. Default = 0.05")
    parser.add_argument("-e", "--entropy", type=str, default = "NSB", help= "Method of entropy estimation. Either NSB or Miller. Default = NSB")
    parser.add_argument('--max', '-x', help="Maximum sample size to calculate the exact entropy of.", type=int)
    parser.add_argument('--logo', help='Produce function logo ps files. If permutation statistical options used the first test in multiple test correction list is used', action="store_true")
    parser.add_argument("-B", help="Number of permutations. Default value is 100", type=int, default=0)
    parser.add_argument("-o", "--stdout", action="store_true", help="Print results to STDOUT")
    parser.add_argument("-M", help = "Specify method to correct p-values for multiple-comparisons. Current methods available: bonferroni, holm, hommel, BH, BY, and hochberg-simes. One or more can be specified separated by colons(:). Default is BY")
    parser.add_argument("struct", help="Structure File")
    parser.add_argument("file_preifx", help="File prefix")
    args = parser.parse_args()

    if (args.text):
        seqStructInfo = seq.SeqStructure(args.struct, "text")
    if (args.cove):
        seqStructInfo = seq.SeqStructure(args.struct, "cove")

    seqStructInfo.parse_sequences(args.file_preifx)

    if (args.max):
        seqStructInfo.calculate_exact(args.max, args.proc)
        if (args.inverse):
            seqStructInfo.calculate_exact(args.max, args.proc, inverse = True)

    if (args.B):
        seqStructInfo.permute(args.B, args.proc)
        print("Calculating permutation information")
        test = seqStructInfo.permInfo(args.entropy, args.proc)
        if (args.inverse):
            print("Calculating inverse permutation information")
            test_inverse = seqStructInfo.permInfo(args.entropy, args.proc, inverse = True)

    if (args.entropy == "NSB"):
        print("Calculating information statistics")
        info, height_dict = seqStructInfo.calculate_entropy_NSB()
        if (args.inverse):
            print("Calculating inverse")
            info_inverse, height_dict_inverse = seqStructInfo.calculate_entropy_inverse_NSB()
    else:
        print("Calculating information statistics using Miller-Maddow estimator")
        info, height_dict = seqStructInfo.calculate_entropy_MM()
        if (args.inverse):
            print("Calculating inverse using MM estimator")
            info_inverse, height_dict_inverse = seqStructInfo.calculate_entropy_inverse_MM()

    if (args.stdout):
        if (args.inverse):
            seqStructInfo.text_output(info, height_dict, info_inverse, height_dict_inverse)
        else:
            seqStructInfo.text_output(info, height_dict)

    if (args.logo):
        seqStructInfo.logo_output(info, height_dict, args.file_preifx, info_inverse, height_dict_inverse)

if __name__ == "__main__":
    main()
