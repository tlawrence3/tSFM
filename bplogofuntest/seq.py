# -*- coding: utf-8 -*-
from collections import Counter, defaultdict
from multiprocessing import Pool
from operator import itemgetter
from string import Template
import bisect
import pkgutil
import itertools
import sys
import numpy as np
import bplogofuntest.nsb_entropy as nb
import random
import time
import math as mt
import bplogofuntest.exact as exact
import glob
import re

class InfoResults:
    def __init__(self):
        infoStates = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))

        self.bpinfodist = defaultdict(int)
        self.bpheightdist = defaultdict(int)

        self.singleinfodist = defaultdict(int)
        self.singleheightdist = defaultdict(int)

    def weighted_dist(self, bpdata, singledata):
        for x in bpdata[0]:
            self.bpinfodist[x] += 1
        for x in bpdata[1]:
            self.bpheightdist[x] += 1

        for x in singledata[0]:
            self.singleinfodist[x] += 1
        for x in singledata[1]:
            self.singleheightdist[x] += 1

    def rtp(data, point, keys_sorted):
        if (point > 0):
            part = 0
            total = sum(data.values())
            i = bisect.bisect_left(keys_sorted, point)
            if (point <= keys_sorted[-1]):
                for y in keys_sorted[i:]:
                    part += data[y]
                return part / total
            else:
                return 0.0
        else:
            return 1.0

class Seq:
    def __init__(self, function, seq):
        self.function = function
        self.seq = seq

    def __len__(self):
        return len(self.seq)


class SeqStructure:
    def __init__(self, struct_file, kind = None, exact = [], inverse = []):
        if (kind):
            self.parse_struct(struct_file, kind)
            self.exact = []
            self.inverse_exact = []
        else:
            self.basepairs = struct_file
            self.exact = exact
            self.inverse_exact = inverse
        self.pos = 0
        self.sequences = []
        self.pairs = set()
        self.singles = set()
        self.functions = Counter()

    def parse_sequences(self, file_prefix):
        for fn in glob.glob("{}_?.aln".format(file_prefix)):
            match = re.search("_([A-Z])\.aln", fn)
            aa_class = match.group(1)
            with open(fn, "r") as ALN:
                good = False
                begin_seq = False
                interleaved = False
                seq = {}
                for line in ALN:
                    match = re.search("^(\S+)\s+(\S+)", line)
                    if (re.search("^CLUSTAL", line)):
                        good = True
                        continue
                    elif (re.search("^[\s\*\.\:]+$", line) and not interleaved and begin_seq):
                        interleaved = True
                    elif (re.search("^[\s\*\.\:]+$", line) and interleaved and begin_seq):
                        continue
                    elif (match and not interleaved):
                        begin_seq = True
                        if (not good):
                            sys.exit("File {} appears not to be a clustal file".format(fn))
                        seq[match.group(1)] = match.group(2)
                    elif (match and interleaved):
                        seq[match.group(1)] += match.group(2)
            for sequence in seq.values():
                self.add_sequence(aa_class, sequence.upper().replace("T", "U"))

        print("{} alignments parsed".format(len(self.functions.keys())), file=sys.stderr)

    def parse_struct(self, struct_file, kind):
        print("Parsing base-pair coordinates", file=sys.stderr)
        basepairs = []
        if (kind == "cove"):
            ss = ""
            pairs = defaultdict(list)
            tarm = 0
            stack = []
            with open(struct_file, "r") as cove:
                for line in cove:
                    line = line.strip()
                    ss += line.split()[1]

            state = "start"
            for count, i in enumerate(ss):
                if (i == ">" and (state == "start" or state == "AD")):
                    if (state == "start"):
                        state = "AD"
                    stack.append(count)

                elif (i == "<" and (state == "AD" or state == "D")):
                    if (state == "AD"):
                        state = "D"
                    pairs[state].append([stack.pop(), count])

                elif (i == ">" and (state == "D" or state == "C")):
                    if (state == "D"):
                        state = "C"
                    stack.append(count)

                elif (i == "<" and (state == "C" or state == "cC")):
                    if (state == "C"):
                        state = "cC"
                    pairs["C"].append([stack.pop(), count])

                elif (i == ">" and (state == "cC" or state == "T")):
                    if (state == "cC"):
                        state = "T"
                    stack.append(count)
                    tarm += 1

                elif (i == "<" and (state == "T" and tarm > 0)):
                    pairs[state].append([stack.pop(), count])
                    tarm -= 1

                elif (i == "<" and (state == "T" or state == "A") and tarm == 0):
                    state = "A"
                    pairs[state].append([stack.pop(), count])

            for arm in pairs:
                for pair in pairs[arm]:
                    basepairs.append((pair[0], pair[1]))


        if (kind == "text"):
             with open(struct_file, "r") as struct:
                for line in struct:
                    line = line.split(";")[1]
                    coords = line.split(",")
                    for coord in coords:
                        pos = coord.split(":")
                        pos = [int(x) for x in pos]
                        basepairs.append((pos[0], pos[1]))

        self.basepairs = basepairs

    def approx_expect(self, H, k, N):
        return H - ((k - 1)/((mt.log(4)) * N))

    def exact_run(self, n, p, numclasses):
        j = exact.calc_exact(n, p, numclasses)
        print("{:2} {:07.5f}".format(n, j[1]), file=sys.stderr)
        return j

    def permuted(self, items, pieces = 2):
        sublists = [[] for i in range(pieces)]
        for x in items:
            sublists[random.randint(0, pieces - 1)].append(x)
        permutedList = []
        for i in range(pieces):
            time.sleep(0.01)
            random.seed()
            random.shuffle(sublists[i])
            permutedList.extend(sublists[i])
        return permutedList

    def permutations(self, numPerm, aa_classes):
        indices = []
        permStructList = []
        print("Generating permuted alignment data", file=sys.stderr)
        for p in range(numPerm):
            indices.append(self.permuted(aa_classes))
        for index in indices:
            permStruct = SeqStructure(self.basepairs, exact = self.exact, inverse = self.inverse_exact)
            for i, seqs in enumerate(self.sequences):
                permStruct.add_sequence(index[i], seqs.seq)
            permStructList.append(permStruct)
        return permStructList

    def permute(self, permute_num, proc):
        with Pool(processes = proc) as pool:
            perm_jobs = []
            for x in range(proc):
                if (x == 0):
                    perm_jobs.append((permute_num//proc+permute_num%proc, self.get_functions()))
                else:
                    perm_jobs.append((permute_num//proc, self.get_functions()))

            perm_results = pool.starmap(self.permutations, perm_jobs)
            self.permutationList = []
            for x in perm_results:
                self.permutationList += x

    def permInfo(self, method, proc, inverse = False):
        with Pool(processes = proc) as pool:
            if (not inverse):
                if (method == "NSB"):
                    perm_info_results = pool.map(self.perm_info_calc_NSB, self.permutationList, len(self.permutationList)//proc)
                else:
                    perm_info_results = pool.map(self.perm_info_calc_MM, self.permutationList, len(self.permutationList)//proc)
            else:
                if (method == "NSB"):
                    perm_info_results = pool.map(self.perm_info_calc_inverse_NSB, self.permutationList, len(self.permutationList)//proc)
                else:
                    perm_info_results = pool.map(self.perm_info_calc_inverse_MM, self.permutationList, len(self.permutationList)//proc)

        return perm_info_results

    def perm_info_calc_MM(self, x):
        total_info_bp = []
        height_info_bp = []
        total_info_ss = []
        height_info_ss = []
        info, height_dict = x.calculate_entropy_MM()

        for coord in sorted(self.basepairs, key = itemgetter(0)):
            if (coord in info):
                for pairtype in sorted(info[coord]):
                    total_info_bp.append(info[coord][pairtype])
                    for aainfo in sorted(height_dict[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                        height_info_bp.append(aainfo[1]*info[coord][pairtype])

        for coord in range(self.pos):
            if (coord in info):
                for base in sorted(info[coord]):
                    total_info_ss.append(info[coord][base])
                    for aainfo in sorted(height_dict[coord][base].items(), key = itemgetter(1), reverse = True):
                        height_info_ss.append(aainfo[1]*info[coord][base])

        return (total_info_bp, total_info_ss, height_info_bp, height_info_ss)

    def perm_info_calc_inverse_MM(self, x):
        total_info_bp = []
        height_info_bp = []
        total_info_ss = []
        height_info_ss = []
        info, height_dict = x.calculate_entropy_inverse_MM()

        for coord in sorted(self.basepairs, key = itemgetter(0)):
            if (coord in info):
                for pairtype in sorted(info[coord]):
                    total_info_bp.append(info[coord][pairtype])
                    for aainfo in sorted(height_dict[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                        height_info_bp.append(aainfo[1]*info[coord][pairtype])

        for coord in range(self.pos):
            if (coord in info):
                for base in sorted(info[coord]):
                    total_info_ss.append(info[coord][base])
                    for aainfo in sorted(height_dict[coord][base].items(), key = itemgetter(1), reverse = True):
                        height_info_ss.append(aainfo[1]*info[coord][base])

        return (total_info_bp, total_info_ss, height_info_bp, height_info_ss)

    def perm_info_calc_inverse_NSB(self, x):
        total_info_bp = []
        height_info_bp = []
        total_info_ss = []
        height_info_ss = []
        info, height_dict = x.calculate_entropy_inverse_NSB()

        for coord in sorted(self.basepairs, key = itemgetter(0)):
            if (coord in info):
                for pairtype in sorted(info[coord]):
                    total_info_bp.append(info[coord][pairtype])
                    for aainfo in sorted(height_dict[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                        height_info_bp.append(aainfo[1]*info[coord][pairtype])

        for coord in range(self.pos):
            if (coord in info):
                for base in sorted(info[coord]):
                    total_info_ss.append(info[coord][base])
                    for aainfo in sorted(height_dict[coord][base].items(), key = itemgetter(1), reverse = True):
                        height_info_ss.append(aainfo[1]*info[coord][base])

        return (total_info_bp, total_info_ss, height_info_bp, height_info_ss)

    def perm_info_calc_NSB(self, x):
        total_info_bp = []
        height_info_bp = []
        total_info_ss = []
        height_info_ss = []
        info, height_dict = x.calculate_entropy_NSB()

        for coord in sorted(self.basepairs, key = itemgetter(0)):
            if (coord in info):
                for pairtype in sorted(info[coord]):
                    total_info_bp.append(info[coord][pairtype])
                    for aainfo in sorted(height_dict[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                        height_info_bp.append(aainfo[1]*info[coord][pairtype])

        for coord in range(self.pos):
            if (coord in info):
                for base in sorted(info[coord]):
                    total_info_ss.append(info[coord][base])
                    for aainfo in sorted(height_dict[coord][base].items(), key = itemgetter(1), reverse = True):
                        height_info_ss.append(aainfo[1]*info[coord][base])

        return (total_info_bp, total_info_ss, height_info_bp, height_info_ss)

    def calculate_exact(self, n, proc, inverse = False):
        exact_list = []
        if (inverse):
            inverse_functions = Counter()
            for aa_class in self.functions:
                inverse_functions[aa_class] = sum(self.functions.values())/self.functions[aa_class]

            p = [x/sum(list(inverse_functions.values())) for x in inverse_functions.values()]
            for i in range(1,n+1):
                exact_list.append((i, p, len(self.functions.values())))

            print("Calculating Sample Size Correction for Inverse", file = sys.stderr)
            with Pool(processes=proc) as pool:
                exact_results = pool.starmap(self.exact_run, exact_list)

            for x in exact_results:
                self.inverse_exact.append(x[1])
        else:
            p = [x/sum(list(self.functions.values())) for x in self.functions.values()]
            for i in range(1,n+1):
                exact_list.append((i, p, len(self.functions.values())))

            print("Calculating Sample Size Correction", file = sys.stderr)
            with Pool(processes=proc) as pool:
                exact_results = pool.starmap(self.exact_run, exact_list)

            for x in exact_results:
                self.exact.append(x[1])

    def calculate_entropy_MM(self):
        info = defaultdict(lambda : defaultdict(float))
        height_dict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))

        functions_array = np.array(list(self.functions.values()))
        bg_entropy = -np.sum((functions_array/functions_array.sum()) * np.log2(functions_array/functions_array.sum()))
        for pairs in self.basepairs:
            for state in self.pairs:
                state_counts = self.get(pairs, state)
                if (sum(state_counts.values()) == 0):
                    continue

                fg_entropy = -np.sum((state_counts/state_counts.sum()) * np.log2(state_counts/state_counts.sum()))
                if (sum(state_counts.values()) <= len(self.exact)):
                    expected_bg_entropy = self.exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = self.approx_expect(bg_entropy, len(self.functions), sum(state_counts.values()))

                if (expected_bg_entropy-fg_entropy < 0):
                    info[pairs][state] = 0
                else:
                    info[pairs][state] = expected_bg_entropy-fg_entropy

                height_class = {}
                for aa_class in state_counts:
                    height_class[aa_class] = (state_counts[aa_class]/sum(state_counts.values())) / (self.functions[aa_class]/len(self))
                for aa_class in height_class:
                    height_dict[pairs][state][aa_class] = height_class[aa_class]/sum(height_class.values())

        for singles in range(self.pos):
            for state in self.singles:
                state_counts = np.array(list(self.get([singles], state).values()))
                if (state_counts.sum() == 0):
                    continue

                fg_entropy = -np.sum((state_counts/state_counts.sum()) * np.log2(state_counts/state_counts.sum()))
                if (sum(state_counts.values()) <= len(self.exact)):
                    expected_bg_entropy = self.exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = self.approx_expect(bg_entropy, len(self.functions), sum(state_counts.values()))

                if (expected_bg_entropy-fg_entropy < 0):
                    info[singles][state] = 0
                else:
                    info[singles][state] = expected_bg_entropy-fg_entropy

                height_class = {}
                for aa_class in state_counts:
                    height_class[aa_class] = (state_counts[aa_class]/sum(state_counts.values())) / (self.functions[aa_class]/len(self))
                for aa_class in height_class:
                    height_dict[singles][state][aa_class] = height_class[aa_class]/sum(height_class.values())

        return (info, height_dict)

    def calculate_entropy_inverse_MM(self):
        info_inverse = defaultdict(lambda : defaultdict(float))
        height_dict_inverse = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
        inverse_functions = Counter()
        for aa_class in self.functions:
            inverse_functions[aa_class] = sum(self.functions.values())/self.functions[aa_class]

        np_inverse_functions = np.array(list(inverse_functions.values()))
        bg_entropy = -np.sum((np_inverse_functions/np_inverse_functions.sum()) * np.log2(np_inverse_functions/np_inverse_functions.sum()))
        for pairs in self.basepairs:
            for state in self.pairs:
                state_counts = self.get(pairs, state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(self.functions)):
                    for function in self.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values())/state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                fg_entropy = -np.sum((nsb_array/nsb_array.sum()) * np.log2(nsb_array/nsb_array.sum()))
                #fg_entropy = nb.S(nb.make_nxkx(nsb_array,nsb_array.size), nsb_array.sum(), nsb_array.size)
                if (sum(state_counts.values()) <= len(self.inverse_exact)):
                    expected_bg_entropy = self.inverse_exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = self.approx_expect(bg_entropy, len(self.functions), sum(state_counts.values()))

                if (expected_bg_entropy-fg_entropy < 0):
                    info_inverse[pairs][state] = 0
                else:
                    info_inverse[pairs][state] = expected_bg_entropy-fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class]/sum(inverse_state_counts.values())) / (inverse_functions[aa_class]/sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[pairs][state][aa_class] = height_class[aa_class]/sum(height_class.values())

        for singles in range(self.pos):
            for state in self.singles:
                state_counts = self.get([singles], state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(self.functions)):
                    for function in self.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values())/state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                fg_entropy = -np.sum((nsb_array/nsb_array.sum()) * np.log2(nsb_array/nsb_array.sum()))
                if (state_counts.sum() <= len(self.inverse_exact)):
                    expected_bg_entropy = self.inverse_exact[state_counts.sum() -1]
                else:
                    expected_bg_entropy = self.approx_expect(bg_entropy, len(self.functions), state_counts.sum())

                if (expected_bg_entropy-fg_entropy < 0):
                    info_inverse[singles][state] = 0
                else:
                    info_inverse[singles][state] = expected_bg_entropy-fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class]/sum(inverse_state_counts.values())) / (inverse_functions[aa_class]/sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[singles][state][aa_class] = height_class[aa_class]/sum(height_class.values())

        return (info_inverse, height_dict_inverse)

    def calculate_entropy_inverse_NSB(self):
        info_inverse = defaultdict(lambda : defaultdict(float))
        height_dict_inverse = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
        inverse_functions = Counter()
        for aa_class in self.functions:
            inverse_functions[aa_class] = sum(self.functions.values())/self.functions[aa_class]

        np_inverse_functions = np.array(list(inverse_functions.values()))
        bg_entropy = -np.sum((np_inverse_functions/np_inverse_functions.sum()) * np.log2(np_inverse_functions/np_inverse_functions.sum()))
        for pairs in self.basepairs:
            for state in self.pairs:
                state_counts = self.get(pairs, state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(self.functions)):
                    for function in self.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values())/state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                fg_entropy = nb.S(nb.make_nxkx(nsb_array,nsb_array.size), nsb_array.sum(), nsb_array.size)
                if (sum(state_counts.values()) <= len(self.inverse_exact)):
                    expected_bg_entropy = self.inverse_exact[sum(state_counts.values()) -1]
                else:
                    expected_bg_entropy = bg_entropy

                if (expected_bg_entropy-fg_entropy < 0):
                    info_inverse[pairs][state] = 0
                else:
                    info_inverse[pairs][state] = expected_bg_entropy-fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class]/sum(inverse_state_counts.values())) / (inverse_functions[aa_class]/sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[pairs][state][aa_class] = height_class[aa_class]/sum(height_class.values())

        for singles in range(self.pos):
            for state in self.singles:
                state_counts = self.get([singles], state)
                if (sum(state_counts.values()) == 0):
                    continue
                if (not len(state_counts) == len(self.functions)):
                    for function in self.functions:
                        state_counts[function] += 1

                inverse_state_counts = Counter()
                for aa_class in state_counts:
                    inverse_state_counts[aa_class] = sum(state_counts.values())/state_counts[aa_class]

                nsb_array = np.array(list(inverse_state_counts.values()))
                fg_entropy = nb.S(nb.make_nxkx(nsb_array,nsb_array.size), nsb_array.sum(), nsb_array.size)
                if (sum(state_counts.values()) <= len(self.inverse_exact)):
                    expected_bg_entropy = self.inverse_exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = bg_entropy

                if (expected_bg_entropy-fg_entropy < 0):
                    info_inverse[singles][state] = 0
                else:
                    info_inverse[singles][state] = expected_bg_entropy-fg_entropy

                height_class = {}
                for aa_class in inverse_state_counts:
                    height_class[aa_class] = (inverse_state_counts[aa_class]/sum(inverse_state_counts.values())) / (inverse_functions[aa_class]/sum(inverse_functions.values()))
                for aa_class in height_class:
                    height_dict_inverse[singles][state][aa_class] = height_class[aa_class]/sum(height_class.values())

        return (info_inverse, height_dict_inverse)

    def calculate_entropy_NSB(self):
        info = defaultdict(lambda : defaultdict(float))
        height_dict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))

        functions_array = np.array(list(self.functions.values()))
        bg_entropy = -np.sum((functions_array/functions_array.sum()) * np.log2(functions_array/functions_array.sum()))
        for pairs in self.basepairs:
            for state in self.pairs:
                state_counts = self.get(pairs, state)
                if (sum(state_counts.values()) == 0):
                    continue
                nsb_array = np.array(list(state_counts.values()) + [0]*(len(self.functions) - len(state_counts)))
                fg_entropy = nb.S(nb.make_nxkx(nsb_array,nsb_array.size), nsb_array.sum(), nsb_array.size)
                if (sum(state_counts.values()) <= len(self.exact)):
                    expected_bg_entropy = self.exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = bg_entropy

                if (expected_bg_entropy-fg_entropy < 0):
                    info[pairs][state] = 0
                else:
                    info[pairs][state] = expected_bg_entropy-fg_entropy

                height_class = {}
                for aa_class in state_counts:
                    height_class[aa_class] = (state_counts[aa_class]/sum(state_counts.values())) / (self.functions[aa_class]/len(self))
                for aa_class in height_class:
                    height_dict[pairs][state][aa_class] = height_class[aa_class]/sum(height_class.values())

        for singles in range(self.pos):
            for state in self.singles:
                state_counts = self.get([singles], state)
                if (sum(state_counts.values()) == 0):
                    continue
                nsb_array = np.array(list(state_counts.values()) + [0]*(len(self.functions) - len(state_counts)))
                fg_entropy = nb.S(nb.make_nxkx(nsb_array,nsb_array.size), nsb_array.sum(), nsb_array.size)
                if (sum(state_counts.values()) <= len(self.exact)):
                    expected_bg_entropy = self.exact[sum(state_counts.values()) - 1]
                else:
                    expected_bg_entropy = bg_entropy

                if (expected_bg_entropy-fg_entropy < 0):
                    info[singles][state] = 0
                else:
                    info[singles][state] = expected_bg_entropy-fg_entropy

                height_class = {}
                for aa_class in state_counts:
                    height_class[aa_class] = (state_counts[aa_class]/sum(state_counts.values())) / (self.functions[aa_class]/len(self))
                for aa_class in height_class:
                    height_dict[singles][state][aa_class] = height_class[aa_class]/sum(height_class.values())

        return (info, height_dict)



    def is_overlap(self, position):
        pass

    def add_sequence(self, function, seq):
        self.sequences.append(Seq(function, seq))
        self.functions[function] += 1
        self.pos = len(seq)
        self.singles.update(seq)
        for x in self.basepairs:
            self.pairs.add(seq[x[0]] +  seq[x[1]])

    def get(self, position, state):
        ret_counter = Counter()
        if (len(position) == 1):
            for x in self.sequences:
                if (x.seq[position[0]] == state[0]):
                    ret_counter[x.function] += 1
        if (len(position) == 2):
            for x in self.sequences:
                if (x.seq[position[0]] == state[0] and x.seq[position[1]] == state[1]):
                    ret_counter[x.function] += 1

        return ret_counter

    def get_functions(self):
        function_list = []
        for key, val in self.functions.items():
            function_list.extend([key]*val)
        return function_list

    def text_output(self, info, height_dict, inverse_info = {}, inverse_height = {}, p = {}):
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
        for coord in sorted(self.basepairs, key = itemgetter(0)):
            if (coord in info):
                for pairtype in sorted(info[coord]):
                    output_string = "bp:\t{}".format(coord)
                    output_string += "\t{}\t{}\t{:05.3f}\t".format(pairtype, sum(self.get(coord, pairtype).values()), info[coord][pairtype])
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

        if (inverse_info):
            print("#ibp\tbp\tN\tinfo{p}{P}".format(**heading_dict))
        for coord in sorted(self.basepairs, key = itemgetter(0)):
            if (coord in inverse_info):
                for pairtype in sorted(inverse_info[coord]):
                    output_string = "ibp:\t{}".format(coord)
                    output_string += "\t{}\t{}\t{:05.3f}\t".format(pairtype, sum(self.get(coord, pairtype).values()), inverse_info[coord][pairtype])
                    if (p):
                        output_string += "{:08.6f}".format(pvalsp[coord][pairtype])
                        for x in multipletesting:
                            output_string += "\t{:08.6f}".format(adjusted_pvals[x]['p'][coord][pairtype])

                    output_string += "\t"
                    for aainfo in sorted(inverse_height[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                        output_string += " {}:{:05.3f}".format(aainfo[0], aainfo[1])
                        if (p):
                            output_string += ":{:08.6f}".format(pvalsP[coord][pairtype][aainfo[0].upper()])
                            for x in multipletesting:
                                output_string += ":{:08.6f}".format(adjusted_pvals[x]['P'][coord][pairtype][aainfo[0].upper()])

                    print(output_string)

        print("#ss\t\tcoord\tf\tN\tinfo{p}{P}".format(**heading_dict))
        for coord in range(self.pos):
            if (coord in info):
                for base in sorted(info[coord]):
                    output_string = "ss:\t\t{}\t{}\t{}\t{:05.3f}".format(coord, base,
                                                                         sum(self.get([coord], base).values()),
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

        if (inverse_info):
            print("#iss\t\tcoord\tf\tN\tinfo{p}{P}".format(**heading_dict))
        for coord in range(self.pos):
            if (coord in inverse_info):
                for base in sorted(inverse_info[coord]):
                    output_string = "iss:\t\t{}\t{}\t{}\t{:05.3f}".format(coord, base,
                                                                         sum(self.get([coord], base).values()),
                                                                         inverse_info[coord][base])
                    if (p):
                        output_string += "\t{:08.6f}".format(pvalsp[coord][base])
                        for x in multipletesting:
                            output_string += "\t{}".format(adjusted_pvals[x]['p'][coord][base])

                    output_string += "\t"
                    for aainfo in sorted(inverse_height[coord][base].items(), key = itemgetter(1), reverse = True):
                        output_string += " {}:{:05.3f}".format(aainfo[0], aainfo[1])
                        if (p):
                            output_string += ":{:08.6f}".format(pvalsP[coord][base][aainfo[0].upper()])
                            for x in multipletesting:
                                output_string += ":{:08.6f}".format(adjusted_pvals[x]['P'][coord][base][aainfo[0].upper()])

                    print(output_string)

    def logo_output(self, info, height_dict, file_prefix, inverse_info = {}, inverse_height = {}, p = {}):
        coord_length = 0 #used to determine eps height
        coord_length_addition = 0

        logo_outputDict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
        inverse_logo_outputDict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))

        #logo output dict construction
        for coord in sorted(self.basepairs, key = itemgetter(0)):
            for pairtype in sorted(self.pairs):
                if (pairtype in info[coord]):
                    for aainfo in sorted(height_dict[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                        logo_outputDict[pairtype][coord][aainfo[0]] = info[coord][pairtype] * aainfo[1]
                else:
                    logo_outputDict[pairtype][coord] = {}

        for coord in range(self.pos):
            for base in sorted(self.singles):
                if (base in info[coord]):
                    for aainfo in sorted(height_dict[coord][base].items(), key = itemgetter(1), reverse = True):
                        logo_outputDict[base][coord][aainfo[0]] = info[coord][base] * aainfo[1]
                else:
                    logo_outputDict[base][coord] = {}

        #inverse logo output dict construction
        for coord in sorted(self.basepairs, key = itemgetter(0)):
            for pairtype in sorted(self.pairs):
                if (pairtype in inverse_info[coord]):
                    for aainfo in sorted(inverse_height[coord][pairtype].items(), key = itemgetter(1), reverse = True):
                        inverse_logo_outputDict[pairtype][coord][aainfo[0]] = inverse_info[coord][pairtype] * aainfo[1]
                else:
                    inverse_logo_outputDict[pairtype][coord] = {}

        for coord in range(self.pos):
            for base in sorted(self.singles):
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
                    if (aainfo[1] < 0.0001 or mt.isnan(aainfo[1])):
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
                    if (aainfo[1] < 0.0001 or mt.isnan(aainfo[1])):
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


    def __len__(self):
        return len(self.sequences)
