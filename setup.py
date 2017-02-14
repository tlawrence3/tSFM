import os
from setuptools import setup, Extension
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
        long_description = f.read()
#bpexact_ext = Extension('bplogofun.exact', ['src/exact.c'])

setup(name = "bplogofuntest",
      install_requires=['statsmodels', 'numpy', 'scipy', 'pandas', 'patsy'],
      packages = ["bplogofuntest"],
      package_data={'bplogofuntest': ['eps/Template.eps', 'eps/inverse_template.eps']},
      entry_points = {
          "console_scripts": ['bplogofuntest = bplogofuntest.bplogofuntest:main']},
      version = "0.1.1",
      description = "Something Something bplogofun",
      long_description=long_description,
      license='GPLv3',
      url = "www.nowhere.com",
      ext_modules=[Extension('bplogofuntest.exact',['src/exact.c'])],)
