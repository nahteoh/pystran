"""
Mechanism for making python aware of the pystran package.
"""

import os
import sys

# When running the tutorials from the command line as
# py .\tutorials\03_weaver_1_tut.py
# the current folder needs to be added to path.
sys.path.insert(0, os.path.abspath("."))
# When running the tutorials from spyder, de tutorial is opened
# in the IDE, and therefore the parent of the current folder needs to be added to path.
sys.path.insert(0, os.path.abspath(".."))
