#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
from collections import defaultdict
from collections import Counter
from operator import itemgetter
from string import Template
from copy import deepcopy
from multiprocessing import Pool
import argparse
import re
import sys
import math as mt
import bplogofuntest.exact
import bplogofuntest.nsb_entropy as nb
import numpy as np
import random
import statsmodels.api
import time
import pkgutil
import bisect
import bplogofuntest.seq as seq

def logo_output(seqStruct, info, height_dict, file_prefix, inverse_info = {}, inverse_height = {}, p = {}):
    coord_length = 0 #used to determine eps height
    coord_length_addition = 0

    logo_outputDict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    inverse_logo_outputDict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))

    #logo output dict construction
    for coord in sorted(seqStruct.basepairs, key = itemgetter(0)):
        for pairtype in sorted(seqStruct.pairs):
            if (pairtype in info[coord]):
                for aainfo in sorted(height_dict[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                    logo_outputDict[pairtype][coord][aainfo[0]] = info[coord][pairtype] * aainfo[1]
            else:
                logo_outputDict[pairtype][coord] = {}

    for coord in range(seqStruct.pos):
        for base in sorted(seqStruct.singles):
            if (base in info[coord]):
                for aainfo in sorted(height_dict[coord][base].items(), key = itemgetter(1), reverse = True):
                    logo_outputDict[base][coord][aainfo[0]] = info[coord][base] * aainfo[1]
            else:
                logo_outputDict[base][coord] = {}

    #inverse logo output dict construction
    for coord in sorted(seqStruct.basepairs, key = itemgetter(0)):
        for pairtype in sorted(seqStruct.pairs):
            if (pairtype in inverse_info[coord]):
                for aainfo in sorted(inverse_height[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                    inverse_logo_outputDict[pairtype][coord][aainfo[0]] = inverse_info[coord][pairtype] * aainfo[1]
            else:
                inverse_logo_outputDict[pairtype][coord] = {}

    for coord in range(seqStruct.pos):
        for base in sorted(seqStruct.singles):
            if (base in inverse_info[coord]):
                for aainfo in sorted(inverse_height[coord][base].items(), key = itemgetter(1), reverse = True):
                    inverse_logo_outputDict[base][coord][aainfo[0]] = inverse_info[coord][base] * aainfo[1]
            else:
                inverse_logo_outputDict[base][coord] = {}

    #output logos
    for base in logo_outputDict:
        logodata = ""
        for coord in sorted(logo_outputDict[base].keys()):
            if (len(str(coord)) > coord_length):
                coord_length = len(str(coord))
            logodata += "numbering {{({}) makenumber}} if\ngsave\n".format(coord)
            for aainfo in sorted(logo_outputDict[base][coord].items(), key = itemgetter(1)):
                if (aainfo[1] < 0.01 or mt.isnan(aainfo[1])):
                    continue
                logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
            logodata += "grestore\nshift\n"
        #output logodata to template
        template_byte = pkgutil.get_data('bplogofuntest', 'eps/Template.eps')
        logo_template = template_byte.decode('utf-8')
        with open("{}_{}.eps".format(base, file_prefix), "w") as logo_output:
            src = Template(logo_template)
            if (len(base) == 2):
                logodata_dict = {'logo_data': logodata, 'low': min(logo_outputDict[base].keys()), 'high': max(logo_outputDict[base].keys()), 'length': 21 * len(logo_outputDict[base].keys()), 'height': 735-(5*(coord_length + coord_length_addition))}
            else:
                logodata_dict = {'logo_data': logodata, 'low': min(logo_outputDict[base].keys()), 'high': max(logo_outputDict[base].keys()), 'length': 15.68 * len(logo_outputDict[base].keys()), 'height': 735-(5*(coord_length + coord_length_addition))}
            logo_output.write(src.substitute(logodata_dict))

    for base in inverse_logo_outputDict:
        logodata = ""
        for coord in sorted(inverse_logo_outputDict[base].keys()):
            if (len(str(coord)) > coord_length):
                coord_length = len(str(coord))
            logodata += "numbering {{({}) makenumber}} if\ngsave\n".format(coord)
            for aainfo in sorted(inverse_logo_outputDict[base][coord].items(), key = itemgetter(1)):
                if (aainfo[1] < 0.01 or mt.isnan(aainfo[1])):
                    continue
                logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
            logodata += "grestore\nshift\n"
        #output logodata to template
        template_byte = pkgutil.get_data('bplogofuntest', 'eps/Template.eps')
        logo_template = template_byte.decode('utf-8')
        with open("inverse_{}_{}.eps".format(base, file_prefix), "w") as logo_output:
            src = Template(logo_template)
            if (len(base) == 2):
                logodata_dict = {'logo_data': logodata, 'low': min(inverse_logo_outputDict[base].keys()), 'high': max(inverse_logo_outputDict[base].keys()), 'length': 21 * len(inverse_logo_outputDict[base].keys()), 'height': 735-(5*(coord_length + coord_length_addition))}
            else:
                logodata_dict = {'logo_data': logodata, 'low': min(inverse_logo_outputDict[base].keys()), 'high': max(inverse_logo_outputDict[base].keys()), 'length': 15.68 * len(inverse_logo_outputDict[base].keys()), 'height': 735-(5*(coord_length + coord_length_addition))}
            logo_output.write(src.substitute(logodata_dict))

def text_output(seqStruct, info, height_dict, inverse_info = {}, inverse_height = {}, p = {}):
    #build output heading
    heading_dict = {}
    if (p):
        pstring = "\tp-value"
        for m in multipletesting:
            pstring += "\t{}".format(m)
        heading_dict['p'] = pstring
    else:
        heading_dict['p'] = ""
    if (p):
        Pstring = "\tclass:height:p-value"
        for m in multipletesting:
            Pstring += ":{}".format(m)
        heading_dict['P'] = Pstring
    else:
        heading_dict['P'] = "\tclass:height"

    print("#bp\tbp\tN\tinfo{p}{P}".format(**heading_dict))
    for coord in sorted(seqStruct.basepairs, key = itemgetter(0)):
        for pairtype in sorted(info[coord]):
            output_string = "bp:\t{}".format(coord)
            output_string += "\t{}\t{}\t{:05.3f}\t".format(pairtype, sum(seqStruct.get(coord, pairtype).values()),
                                                       info[coord][pairtype])
            if (p):
                output_string += "{:08.6f}".format(pvalsp[coord][pairtype])
                for x in multipletesting:
                    output_string += "\t{:08.6f}".format(adjusted_pvals[x]['p'][coord][pairtype])

            output_string += "\t"
            for aainfo in sorted(height_dict[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                output_string += " {}:{:05.3f}".format(aainfo[0], aainfo[1])
                if (p):
                    output_string += ":{:08.6f}".format(pvalsP[coord][pairtype][aainfo[0].upper()])
                    for x in multipletesting:
                        output_string += ":{:08.6f}".format(adjusted_pvals[x]['P'][coord][pairtype][aainfo[0].upper()])

            print(output_string)

    print("#ss\t\tcoord\tf\tN\tinfo{p}{P}".format(**heading_dict))
    for coord in range(seqStruct.pos):
        for base in sorted(info[coord]):
            output_string = "ss:\t\t{}\t{}\t{}\t{:05.3f}".format(coord, base,
                                                                 sum(seqStruct.get([coord], base).values()),
                                                                 info[coord][base])
            if (p):
                output_string += "\t{:08.6f}".format(pvalsp[coord][base])
                for x in multipletesting:
                    output_string += "\t{}".format(adjusted_pvals[x]['p'][coord][base])

            output_string += "\t"
            for aainfo in sorted(height_dict[coord][base].items(), key = itemgetter(1), reverse = True):
                output_string += " {}:{:05.3f}".format(aainfo[0], aainfo[1])
                if (p):
                    output_string += ":{:08.6f}".format(pvalsP[coord][base][aainfo[0].upper()])
                    for x in multipletesting:
                        output_string += ":{:08.6f}".format(adjusted_pvals[x]['P'][coord][base][aainfo[0].upper()])

            print(output_string)

def main():
     #Setup parser
    parser = argparse.ArgumentParser(description = "bpLogoFun")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--infernal", action="store_true", help="Structure file is in infernal format")
    group.add_argument("-c", "--cove", action="store_true", help="Structure file is in cove format")
    group.add_argument("-t", "--text", action="store_true", help="Structure file is in text format")
    parser.add_argument("-v","--inverse", action="store_true", help="calculate anti-determinates")
    parser.add_argument("-a", "--alpha", type=float, default=0.05, help="Alpha value used for statistical tests. Default = 0.05")
    parser.add_argument('--max', '-x', help="Maximum sample size to calculate the exact entropy of.", type=int)
    parser.add_argument('--logo', help='Produce function logo ps files. If permutation statistical options used the first test in multiple test correction list is used', action="store_true")
    parser.add_argument("-s", "--single", action="store_true", help="Calculate information statistics for single sites in addition to base-pair information statistics.")
    parser.add_argument("-B", help="Number of permutations. Default value is 100", type=int, default=0)
    parser.add_argument("-o", "--stdout", action="store_true", help="Print results to STDOUT")
    parser.add_argument("-M", default="BY", help = "Specify method to correct p-values for multiple-comparisons. Current methods available: bonferroni, holm, hommel, BH, BY, and hochberg-simes. One or more can be specified separated by colons(:). Default is BY")
    parser.add_argument("-d", action="store_true", help="Output the alignment coordinates that correspond to each base-pair")
    parser.add_argument("struct", help="Structure File")
    parser.add_argument("file_preifx", help="File prefix")
    args = parser.parse_args()

    if (args.text):
        seqStructInfo = seq.SeqStructure(args.struct, "text")
    if (args.cove):
        seqStructInfo = seq.SeqStructure(args.struct, "cove")

    seqStructInfo.parse_sequences(args.file_preifx)

    if (args.B):
        seqStructInfo.permute(args.B)

    if (args.max):
        seqStructInfo.calculate_exact(args.x)

    info = defaultdict(lambda : defaultdict(float))
    height_dict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))

    print("Calculating information statistics")
    functions_array = np.array(list(seqStructInfo.functions.values()))
    bg_entropy = -np.sum((functions_array/functions_array.sum()) * np.log2(functions_array/functions_array.sum()))
    for pairs in seqStructInfo.basepairs:
        for state in seqStructInfo.pairs:
            state_counts = seqStructInfo.get(pairs, state)
            if (sum(state_counts.values()) == 0):
                continue
            nsb_array = np.array(list(state_counts.values()) + [0]*(len(seqStructInfo.functions) - len(state_counts)))
            fg_entropy = nb.S(nb.make_nxkx(nsb_array,nsb_array.size), nsb_array.sum(), nsb_array.size)
            if (bg_entropy-fg_entropy < 0):
                info[pairs][state] = 0
            else:
                info[pairs][state] = bg_entropy-fg_entropy

            height_class = {}
            for aa_class in state_counts:
                height_class[aa_class] = (state_counts[aa_class]/sum(state_counts.values())) / (seqStructInfo.functions[aa_class]/len(seqStructInfo))
            for aa_class in height_class:
                height_dict[pairs][state][aa_class] = height_class[aa_class]/sum(height_class.values())

    for singles in range(seqStructInfo.pos):
        for state in seqStructInfo.singles:
            state_counts = seqStructInfo.get([singles], state)
            if (sum(state_counts.values()) == 0):
                continue
            nsb_array = np.array(list(state_counts.values()) + [0]*(len(seqStructInfo.functions) - len(state_counts)))
            fg_entropy = nb.S(nb.make_nxkx(nsb_array,nsb_array.size), nsb_array.sum(), nsb_array.size)
            if (bg_entropy-fg_entropy < 0):
                info[singles][state] = 0
            else:
                info[singles][state] = bg_entropy-fg_entropy

            height_class = {}
            for aa_class in state_counts:
                height_class[aa_class] = (state_counts[aa_class]/sum(state_counts.values())) / (seqStructInfo.functions[aa_class]/len(seqStructInfo))
            for aa_class in height_class:
                height_dict[singles][state][aa_class] = height_class[aa_class]/sum(height_class.values())

    if (args.inverse):
        print("Calculating inverse")
        info_inverse = defaultdict(lambda : defaultdict(float))
        height_dict_inverse = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
        inverse_functions = Counter()
        for aa_class in seqStructInfo.functions:
            inverse_functions[aa_class] = sum(seqStructInfo.functions.values())/seqStructInfo.functions[aa_class]

        np_inverse_functions = np.array(list(inverse_functions.values()))
        bg_entropy = -np.sum((np_inverse_functions/np_inverse_functions.sum()) * np.log2(np_inverse_functions/np_inverse_functions.sum()))
        for pairs in seqStructInfo.basepairs:
            for state in seqStructInfo.pairs:
                state_counts = seqStructInfo.get(pairs, state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(seqStructInfo.functions)):
                    for function in seqStructInfo.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values())/state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                fg_entropy = nb.S(nb.make_nxkx(nsb_array,nsb_array.size), nsb_array.sum(), nsb_array.size)
                if (bg_entropy-fg_entropy < 0):
                    info_inverse[pairs][state] = 0
                else:
                    info_inverse[pairs][state] = bg_entropy-fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class]/sum(inverse_state_counts.values())) / (inverse_functions[aa_class]/sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[pairs][state][aa_class] = height_class[aa_class]/sum(height_class.values())

        for singles in range(seqStructInfo.pos):
            for state in seqStructInfo.singles:
                state_counts = seqStructInfo.get([singles], state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(seqStructInfo.functions)):
                    for function in seqStructInfo.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values())/state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                fg_entropy = nb.S(nb.make_nxkx(nsb_array,nsb_array.size), nsb_array.sum(), nsb_array.size)
                if (bg_entropy-fg_entropy < 0):
                    info_inverse[singles][state] = 0
                else:
                    info_inverse[singles][state] = bg_entropy-fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class]/sum(inverse_state_counts.values())) / (inverse_functions[aa_class]/sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[singles][state][aa_class] = height_class[aa_class]/sum(height_class.values())

    if (args.stdout):
        text_output(seqStructInfo, info, height_dict)

    if (args.logo):
        logo_output(seqStructInfo, info, height_dict, args.file_preifx, info_inverse, height_dict_inverse)

if __name__ == "__main__":
    main()
