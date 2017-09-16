import tsfm.MolecularInformation

def test_functionlogo_cove(cove_files):
    basepair_list = [(12, 24), (11, 25), (10, 26), (9, 27), (33, 41), (32, 42), (31, 43),
                    (30, 44), (29, 45), (55, 63), (54, 64), (53, 65), (52, 66), (51, 67),
                    (6, 68), (5, 69), (4, 70), (3, 71), (2, 72), (1, 73), (0, 74)]
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    assert sorted(basepair_list) == sorted(cove_logo.basepairs)

def test_parseSequence(cove_files):
    cove_logo = tsfm.MolecularInformation.FunctionLogo(cove_files['cove'], "cove")
    cove_logo.parse_sequences(cove_files['prefix'])
    assert set(['UA', 'CC', 'CA', 'AA', 'CG', 'GC', 'AU', 'GU']) == cove_logo.pairs
