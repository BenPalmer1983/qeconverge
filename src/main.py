################################################################
#    Processing PWscf input file
#
#
#
#
################################################################

import os
import datetime
import re
import sys
from pwscf_converge import pwscf_converge


new_run = pwscf_converge()
new_run.run()
