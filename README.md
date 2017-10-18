[![Travis](https://img.shields.io/travis/tlawrence3/tsfm.svg)]() [![PyPI](https://img.shields.io/pypi/format/tsfm.svg)]() [![PyPI](https://img.shields.io/pypi/pyversions/tsfm.svg)]() [![PyPI](https://img.shields.io/pypi/l/tsfm.svg)]()
[![Coverage Status](https://coveralls.io/repos/github/tlawrence3/tsfm/badge.svg?branch=master)](https://coveralls.io/github/tlawrence3/tsfm?branch=master)
# Read the docs
http://tsfm-trna-structure-function-mapper.readthedocs.io/en/latest/
# Example command
```shell
tsfm -vc struct_file.txt ferns/FERN Gnetidae/GNET -j
```
This example command assumes you have the below directory structure and follow this file naming convention (`otuName_function.aln`) for sequence alignment files in clustal format:
```shell
.
├── Gnetidae
│   ├── GNET_A.aln
│   ├── GNET_C.aln
│   ├── GNET_D.aln
│   ├── GNET_E.aln
│   ├── GNET_F.aln
│   ├── GNET_G.aln
│   ├── GNET_H.aln
│   ├── GNET_I.aln
│   ├── GNET_J.aln
│   ├── GNET_K.aln
│   ├── GNET_L.aln
│   ├── GNET_M.aln
│   ├── GNET_N.aln
│   ├── GNET_P.aln
│   ├── GNET_Q.aln
│   ├── GNET_R.aln
│   ├── GNET_S.aln
│   ├── GNET_T.aln
│   ├── GNET_V.aln
│   ├── GNET_W.aln
│   ├── GNET_X.aln
│   └── GNET_Y.aln
└── ferns
    ├── FERN_A.aln
    ├── FERN_C.aln
    ├── FERN_D.aln
    ├── FERN_E.aln
    ├── FERN_F.aln
    ├── FERN_G.aln
    ├── FERN_H.aln
    ├── FERN_I.aln
    ├── FERN_J.aln
    ├── FERN_K.aln
    ├── FERN_L.aln
    ├── FERN_M.aln
    ├── FERN_N.aln
    ├── FERN_P.aln
    ├── FERN_Q.aln
    ├── FERN_R.aln
    ├── FERN_S.aln
    ├── FERN_T.aln
    ├── FERN_V.aln
    ├── FERN_W.aln
    ├── FERN_X.aln
    └── FERN_Y.aln
```
