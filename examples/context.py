import os
import sys

# This is an ugly hack. Jupyter notebook runs in a different path than plain python.
sys.path.insert(0, os.path.abspath('../stranalyzer'))
sys.path.insert(0, os.path.abspath('../stranalyzer/stranalyzer'))

import stranalyzer