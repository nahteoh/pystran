"""
Mechanism for making python aware of the pystran package.
"""

import os
import sys

# This is an ugly hack. Jupyter notebook runs in a different path than plain python.
sys.path.insert(0, os.path.abspath("../.."))
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("."))

import pystran
