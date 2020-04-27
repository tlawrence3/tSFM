# -*- coding: utf-8 -*-
"""
This module contains tests for the MolecularInformation module.
"""

from pytest import approx
import tsfm.MolecularInformation

def test_functionlogo_text_setup(cove_files):
    """
    Testing initialization of function logo object using text format
    for basepair annotation.
    """
    basepair_list = [(0, 72), (1, 71), (2, 70), (3, 69), (4, 68), (5, 67),
                     (6, 66), (9, 25), (10, 24), (11, 23), (12, 22), (27, 43),
                     (28, 42), (29, 41), (30, 40), (31, 39), (49, 65), (50, 64),
                     (51, 63), (52, 62), (53, 61)]
    text_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['text'], "text")
    assert sorted(basepair_list) == sorted(text_logo.basepairs)

def test_functionlogo_cove_setup(cove_files):
    """
    Testing initialization of function logo object using coves format
    for basepair annotation.
    """
    basepair_list = [(12, 24), (11, 25), (10, 26), (9, 27), (33, 41), (32, 42), (31, 43),
                     (30, 44), (29, 45), (55, 63), (54, 64), (53, 65), (52, 66), (51, 67),
                     (6, 68), (5, 69), (4, 70), (3, 71), (2, 72), (1, 73), (0, 74)]
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    assert sorted(basepair_list) == sorted(cove_logo.basepairs)
    assert set(['UA', 'CC', 'CA', 'AA', 'CG', 'GC', 'AU', 'GU']) == cove_logo.pairs
    assert len(cove_logo) == 21
    assert cove_logo.pos == 76
    assert set(['-', 'A', 'U', 'C', 'G']) == cove_logo.singles
    assert {'H': 14, 'K': 7} == cove_logo.functions

def test_functionlogo_exact1(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.calculate_exact(5, 1)
    cove_logo.calculate_exact(5, 1, inverse=True)
    assert len(cove_logo.exact) == 5
    assert approx([0.0, 0.44442903, 0.61221003, 0.69694880, 0.746884728]) == cove_logo.exact
    assert len(cove_logo.inverse_exact) == 5
    assert approx([0.0, 0.444429036427, 0.612210037675, 0.696948806402, 0.746884728233]) == cove_logo.inverse_exact

def test_functionlogo_exact2(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.calculate_exact(5, 10)
    cove_logo.calculate_exact(5, 10, inverse=True)
    assert len(cove_logo.exact) == 5
    assert approx([0.0, 0.44442903, 0.61221003, 0.69694880, 0.746884728]) == cove_logo.exact
    assert len(cove_logo.inverse_exact) == 5
    assert approx([0.0, 0.444429036427, 0.612210037675, 0.696948806402, 0.746884728233]) == cove_logo.inverse_exact

def test_functionlogo_cove_MM(cove_files):
    info_key1 = [(12, 24), (11, 25), (10, 26), (9, 27), (33, 41), (32, 42),
                 (31, 43), (30, 44), (29, 45), (55, 63), (54, 64), (53, 65),
                 (52, 66), (51, 67), (6, 68), (5, 69), (4, 70), (3, 71),
                 (2, 72), (1, 73), (0, 74), 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
                 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
                 71, 72, 73, 74, 75]
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.calculate_exact(5, 1)
    cove_logo.calculate_exact(5, 1, inverse=True)
    info, height = cove_logo.calculate_entropy_MM()
    inverse_info, inverse_height = cove_logo.calculate_entropy_inverse_MM()
    assert sorted(info_key1, key=str) == sorted(list(info.keys()), key=str)
    assert approx({'U': 0.8152461882767065, 'C': 0.866771011165598}) == info[27]
    assert sorted(info_key1, key=str) == sorted(list(height.keys()), key=str)
    assert sorted(['UA', 'GU']) == sorted(list(height[(3, 71)].keys()))
    assert approx({'K': 1.0}) == height[(3, 71)]['UA']
    assert sorted(info_key1, key=str) == sorted(list(inverse_info.keys()), key=str)
    assert sorted(info_key1, key=str) == sorted(list(inverse_height.keys()), key=str)
    assert approx({'UA': 0.33488777478501253, 'CG': 0.53592154740969555}) == inverse_info[(11, 25)]
    assert approx({'H': 0.11764705882352941, 'K': 0.8823529411764706}) == inverse_height[(11, 25)]['CG']

def test_functionlogo_cove_NSB(cove_files):
    info_key1 = [(12, 24), (11, 25), (10, 26), (9, 27), (33, 41), (32, 42),
                 (31, 43), (30, 44), (29, 45), (55, 63), (54, 64), (53, 65),
                 (52, 66), (51, 67), (6, 68), (5, 69), (4, 70), (3, 71),
                 (2, 72), (1, 73), (0, 74), 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
                 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
                 71, 72, 73, 74, 75]
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.calculate_exact(5, 1)
    cove_logo.calculate_exact(5, 1, inverse=True)
    info, height = cove_logo.calculate_entropy_NSB()
    inverse_info, inverse_height = cove_logo.calculate_entropy_inverse_NSB()
    assert sorted(info_key1, key=str) == sorted(list(info.keys()), key=str)
    assert sorted(info_key1, key=str) == sorted(list(height.keys()), key=str)
    assert sorted(info_key1, key=str) == sorted(list(inverse_info.keys()), key=str)
    assert sorted(info_key1, key=str) == sorted(list(inverse_height.keys()), key=str)
    assert approx({'A': 0.0054504094063251296}) == info[60]
    assert approx({'G': 0.046427861927731184}) == inverse_info[32]
    assert sorted(['CG', 'GC']) == sorted(list(height[(6, 68)].keys()))
    assert sorted(['CG', 'GC']) == sorted(list(inverse_height[(6, 68)].keys()))
    assert approx({'K': 1.0}) == height[(6, 68)]['CG']
    assert approx({'K': 0.05882352941176471, 'H': 0.9411764705882354}) == inverse_height[(6, 68)]['CG']

def test_functionlogo_cove_perm1(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.permute(4, 1)
    perm_info_miller = cove_logo.permInfo("Miller", 1)
    perm_info_nsb = cove_logo.permInfo("NSB", 1)
    perm_inverse_miller = cove_logo.permInfo("Miller", 1, inverse=True)
    perm_inverse_nsb = cove_logo.permInfo("NSB", 1, inverse=True)

def test_functionlogo_cove_perm2(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.permute(4, 8)
    perm_info_miller = cove_logo.permInfo("Miller", 8)
    perm_info_nsb = cove_logo.permInfo("NSB", 8)
    perm_inverse_miller = cove_logo.permInfo("Miller", 8, inverse=True)
    perm_inverse_nsb = cove_logo.permInfo("NSB", 8, inverse=True)

def test_FunctionLogoResult(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    result = tsfm.MolecularInformation.FunctionLogoResults(cove_files['prefix'],
                                                           cove_logo.basepairs,
                                                           cove_logo.pos,
                                                           cove_logo.sequences,
                                                           cove_logo.pairs,
                                                           cove_logo.singles)

def test_statTest(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.permute(4, 8)
    perm_info_miller = cove_logo.permInfo("Miller", 8)
    perm_info_nsb = cove_logo.permInfo("NSB", 8)
    perm_inverse_miller = cove_logo.permInfo("Miller", 8, inverse=True)
    perm_inverse_nsb = cove_logo.permInfo("NSB", 8, inverse=True)
    result_miller = tsfm.MolecularInformation.FunctionLogoResults(cove_files['prefix'],
                                                                  cove_logo.basepairs,
                                                                  cove_logo.pos,
                                                                  cove_logo.sequences,
                                                                  cove_logo.pairs,
                                                                  cove_logo.singles)
    result_nsb = tsfm.MolecularInformation.FunctionLogoResults(cove_files['prefix'],
                                                               cove_logo.basepairs,
                                                               cove_logo.pos,
                                                               cove_logo.sequences,
                                                               cove_logo.pairs,
                                                               cove_logo.singles)
    info_nsb, height_nsb = cove_logo.calculate_entropy_NSB()
    inverse_info_nsb, inverse_height_nsb = cove_logo.calculate_entropy_inverse_NSB()
    info_miller, height_miller = cove_logo.calculate_entropy_MM()
    inverse_info_miller, inverse_height_miller = cove_logo.calculate_entropy_inverse_MM()

    result_miller.add_information(info=info_miller, height=height_miller)
    result_miller.add_information(info=inverse_info_miller, height=inverse_height_miller, inverse=True)
    result_nsb.add_information(info=info_nsb, height=height_nsb)
    result_nsb.add_information(info=inverse_info_nsb, height=inverse_height_nsb, inverse=True)
    result_miller.add_stats(perm_info_miller, 'fdr_bh')
    result_miller.add_stats(perm_inverse_miller, 'fdr_bh', inverse=True)
    result_nsb.add_stats(perm_info_nsb, 'fdr_bh')
    result_nsb.add_stats(perm_inverse_nsb, 'fdr_bh', inverse=True)
