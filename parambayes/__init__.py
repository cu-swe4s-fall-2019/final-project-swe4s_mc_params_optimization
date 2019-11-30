"""
ParamBayes
Bayesian MCMC parameterization for CSCI 7000 (swe4s) class project
"""

# Add imports here
print(__file__)
from .parambayes import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
