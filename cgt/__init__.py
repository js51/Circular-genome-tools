"""
"""
name = "cgt"
from sage.all_cmdline import *
gigs = 30
old_command = gap._Expect__command
s = old_command.index('-o ')
e = old_command.index(' -', s)
gap._Expect__command = gap._Expect__command.replace(gap._Expect__command[s:e], f'-o {gigs}G')
from . import visualisation
from . import models
from . import position_paradigm
from . import enums
from . import rearrangements
from . import examples
from .examples import *
from .enums import *
from .position_paradigm import PositionParadigmFramework
from .models import Model