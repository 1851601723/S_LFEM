# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 16:50:42 2021

@author: dell
"""

from os import getcwd
import sys
path = getcwd()
path = list(path)
path = path + ['\\','l','c','m']
path = ''.join(path)
sys.path.append(path)

path = getcwd()
path = list(path)
path = path + ['\\','i','n','p']
path = ''.join(path)
sys.path.append(path)

import method
import inpNeoHookeanBeam as inp1

a = method.mNewtonStatics(inp1)

