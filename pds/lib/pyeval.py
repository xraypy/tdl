"""
Simple python interpretor functions

Authors / Modifications:
------------------------
* T. Trainor (tptrainor@alaska.edu) 10-6-2008  
* Original version: numeval.py from data_shell_2 (6-1-2004)

Description:
------------ 
This module just sets up a clean name-space for execution
of python commands/code.  VAR_DICT should be a dictionary
defined in another module, that keep track of it.  E.g. Shell

Note that exception handling is not done in the below exec/eval
calls intentionally.  If there is an exception it should be raised
so the caller knows about it

Note that imports are limited to keep a clean
namespace for eval functions.  ie limit to
minimum the code in this module....
 
"""
###################################################################

import types

###################################################################
def init_namespace():
    """
    Create a variables dictionary for code exec
    delete stuff not needed....
    """
    vars = {}
    tmp  = globals()
    vars = tmp.copy()
    vars["__name__"] = "__main__"
    del vars['__doc__']
    del vars['__file__']
    del vars['init_namespace']
    del vars['do_exec']
    del vars['do_execfile']
    del vars['do_eval']
    del vars['types']
    return vars

###################################################################
def do_exec(exec_str,vars={}):
    """
    Exec code str in vars workspace
    return 0 if command is completed

    Note this function does not deal with
    line continuations.  Use code.InteractiveConsole
    module for full python emulation
    """
    if vars == {}: vars = init_namespace()
    if type(exec_str) != types.ListType: exec_str = [exec_str]
    #exec(exec_str,vars)  # ,vars)
    for s in exec_str:
        exec s in vars
    return 0

###################################################################
def do_execfile(fname,vars={}):
    """
    Exec file in vars workspace 
    """
    if vars == {}: vars = init_namespace()
    execfile(fname, vars)
    return 0

###################################################################
def do_eval(eval_str,vars={}):
    """
    Do an eval and ret data
    """
    if vars == {}: vars = init_namespace()
    ret = eval(eval_str,vars)
    return ret

###################################################################
###################################################################
# Test
if __name__ == "__main__":
    v = init_namespace()
    print v, "\n----\n"
    
    do_exec("x=10;y=100",vars=v)
    print v, "\n----\n"

    s = "def t():\n   return 100\nxx = t()"
    do_exec(s,vars=v)
    s = "import numpy as num"
    do_exec(s,vars=v)
    print v, "\n----\n"
    
