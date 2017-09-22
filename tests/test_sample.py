import tsfm.MolecularInformation
from pytest import approx

def test_functionlogo_text_setup(cove_files):
    basepair_list =[(0, 72), (1, 71), (2, 70), (3, 69), (4, 68), (5, 67),
                    (6, 66), (9, 25), (10, 24), (11, 23), (12, 22), (27, 43),
                    (28, 42), (29, 41), (30, 40), (31, 39), (49, 65), (50, 64),
                    (51, 63), (52, 62), (53, 61)]
    text_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['text'], "text")
    assert sorted(basepair_list) == sorted(text_logo.basepairs) 

def test_functionlogo_cove_setup(cove_files):
    basepair_list = [(12, 24), (11, 25), (10, 26), (9, 27), (33, 41), (32, 42), (31, 43),
                    (30, 44), (29, 45), (55, 63), (54, 64), (53, 65), (52, 66), (51, 67),
                    (6, 68), (5, 69), (4, 70), (3, 71), (2, 72), (1, 73), (0, 74)]
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    assert sorted(basepair_list) == sorted(cove_logo.basepairs)
    assert set(['UA', 'CC', 'CA', 'AA', 'CG', 'GC', 'AU', 'GU']) == cove_logo.pairs
    assert 21 == len(cove_logo)
    assert 76 == cove_logo.pos
    assert set(['-', 'A', 'U', 'C', 'G']) == cove_logo.singles
    assert {'H': 14, 'K': 7} == cove_logo.functions

def test_functionlogo_exact1(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.calculate_exact(5, 1)
    cove_logo.calculate_exact(5,1,inverse=True)
    assert len(cove_logo.exact) == 5
    assert approx([0.0, 0.44442903, 0.61221003, 0.69694880, 0.746884728]) == cove_logo.exact
    assert len(cove_logo.inverse_exact) == 5
    assert approx([0.0, 0.444429036427, 0.612210037675, 0.696948806402, 0.746884728233]) == cove_logo.inverse_exact

def test_functionlogo_exact2(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.calculate_exact(5, 10)
    cove_logo.calculate_exact(5,10,inverse=True)
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
    cove_logo.calculate_exact(5,1,inverse=True)
    info, height = cove_logo.calculate_entropy_MM()
    inverse_info, inverse_height = cove_logo.calculate_entropy_inverse_MM()
    assert info_key1 == list(info.keys())
    assert approx(sorted({'U': 0.8152461882767065, 'C': 0.866771011165598})) == sorted(info[27])
    assert info_key1 == list(height.keys())
    assert sorted(['UA', 'GU']) == sorted(list(height[(3, 71)].keys()))
    assert approx({'K': 1.0}) == height[(3, 71)]['UA']
    assert info_key1 == list(inverse_info.keys())
    assert info_key1 == list(inverse_height.keys())
    assert approx(sorted({'UA': 0.33488777478501253, 'CG': 0.53592154740969555})) == sorted(inverse_info[(11, 25)])
    assert approx(sorted({'H': 0.11764705882352941, 'K': 0.8823529411764706})) == sorted(inverse_height[(11, 25)]['CG'])

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
    cove_logo.calculate_exact(5,1,inverse=True)
    info, height = cove_logo.calculate_entropy_NSB()
    inverse_info, inverse_height = cove_logo.calculate_entropy_inverse_NSB()
    assert info_key1 == list(info.keys())
    assert info_key1 == list(height.keys())
    assert info_key1 == list(inverse_info.keys())
    assert info_key1 == list(inverse_height.keys())
    assert approx({'A': 0.0054504094063251296}) == info[60]
    assert approx({'G': 0.046427861927731184}) == inverse_info[32]
    assert sorted(['CG', 'GC']) == sorted(list(height[(6, 68)].keys()))
    assert sorted(['CG', 'GC']) == sorted(list(inverse_height[(6, 68)].keys()))
    assert approx(sorted({'K': 1.0})) == sorted(height[(6, 68)]['CG'])
    assert approx(sorted({'K': 0.05882352941176471, 'H': 0.9411764705882354})) == sorted(inverse_height[(6, 68)]['CG'])

def test_functionlogo_cove_perm1(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.permute(50, 1)

def test_functionlogo_cove_perm2(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.permute(5, 8)
