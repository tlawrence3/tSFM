import tsfm.MolecularInformation
from pytest import approx

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

def test_functionlogo_cove_exact1(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.calculate_exact(5, 1)
    cove_logo.calculate_exact(5,1,inverse=True)
    assert len(cove_logo.exact) == 5
    assert approx([0.0, 0.44442903, 0.61221003, 0.69694880, 0.746884728]) == cove_logo.exact
    assert len(cove_logo.inverse_exact) == 5
    assert approx([0.0, 0.444429036427, 0.612210037675, 0.696948806402, 0.746884728233]) == cove_logo.inverse_exact

def test_functionlogo_cove_exact2(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.calculate_exact(5, 10)
    cove_logo.calculate_exact(5,10,inverse=True)
    assert len(cove_logo.exact) == 5
    assert approx([0.0, 0.44442903, 0.61221003, 0.69694880, 0.746884728]) == cove_logo.exact
    assert len(cove_logo.inverse_exact) == 5
    assert approx([0.0, 0.444429036427, 0.612210037675, 0.696948806402, 0.746884728233]) == cove_logo.inverse_exact

def test_functionlogo_cove_MM(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    cove_logo.calculate_exact(5, 1)
    cove_logo.calculate_exact(5,1,inverse=True)
    info, height = cove_logo.calculate_entropy_MM()
