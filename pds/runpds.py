"""
Start the pds shell program

Notes:
------
Exectute this script as:
>>python runpds.py
- or- 
>>python runpds.py args

To see the startup options:
>>python runpds.py -h

This script can also be imported at the python prompt
(if its on your path):
>>python
>>import runpds
>>runpds.main()

"""
def main(args=''):
    """
    Start the shell

    This first checks to see if pds is on the path.
    If not we add the root directory so 'import pds'
    will work.  
    """
    try:
        from pds import shell
        shell.main(args)
    except:
        import os, sys
        root = os.path.dirname(os.path.abspath(__file__))
        root = os.path.dirname(root)
        sys.path.append(root)
        #print '--', __file__,',', root
        #for pp in sys.path: print '#',pp
        from pds import shell
        shell.main(args)
    
#####################################################################################
if (__name__ == '__main__'):
    import sys, os
    args = sys.argv[1:]
    if len(args) > 0:
        if args[0].lower() == 'pds':
            args.pop(0)
    #print args
    main(args)
