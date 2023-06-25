"""
"""
name = "cgt"
from sage.all_cmdline import gap
gigs = 30
old_command = gap._Expect__command
try:
    s = old_command.index('-o ')
    e = old_command.index(' -', s)
    gap._Expect__command = gap._Expect__command.replace(gap._Expect__command[s:e], f'-o {gigs}G')
except ValueError:
    s = old_command.index('-b ')
    e = old_command.index(' -', s)
    gap._Expect__command = gap._Expect__command.replace(gap._Expect__command[s:e], f'-b -o {gigs}G')
from .examples import *
from .enums import *
from .position_paradigm import PositionParadigmFramework
from .models import Model