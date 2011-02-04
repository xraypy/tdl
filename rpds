#!python
"""
Startup pds
"""
##############################################################
import sys
args = sys.argv[1:]
pds = True
if len(args) > 0:
    if args[0].lower() == 'pds':
        pds = True
        args.pop(0)

if pds:
    from tdl import pds
    pds.main(args)
else:
    pass
