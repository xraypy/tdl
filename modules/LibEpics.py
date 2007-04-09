#!/usr/bin/env python
##
## tdl interface to EpicsCA for Epics Control system

import EpicsCA
v = EpicsCA.__version__

class _tdl_cb:
    def __init__(self,pv=None,cb=None,tdl=None):
        self.kws = {'pv':pv}
        self.cb  = cb
        self.tdl = tdl
    def __call__(self,**kws):
        self.kws.update(kws)
        self.tdl.run_procedure(self.cb,args=self.cb.args,kws=self.kws)

def _set_callback(pv=None,callback=None,tdl=None,**kw):
    return pv.add_callback(callback=_tdl_cb(pv=pv, tdl=tdl,
                                            cb=tdl.getVariable(callback)))

def _pv(pvname, **kw):
    """connect to an Epics PV:
       pv = _epics.pv('XXX:m1.VAL')
       print pv.pvname, pv.value
       pv.put(1)
    """
    return EpicsCA.PV(pvname, **kw)

def _clear_callback(pv):              return pv.clear_callbacks()
def _caput(pv,val,**kw):              return EpicsCA.caput(pv,val,**kw)
def _caget(x,use_char=False, **kw):   return EpicsCA.caget(x,use_char=use_char,**kw)
def _pend_event(t=1.e-2):             return EpicsCA.pend_event(t)
def _pend_io(t=1.e-2):                return EpicsCA.pend_io(t)

HelpEpics = """
    Epics PV:  simple interface to Epics Channel Access

"""

_groups_ = [('_epics',True)]
_help_   = {'epics': HelpEpics}
_func_   = {'_epics.pv':(_pv, None),
            '_epics.set_callback':(_set_callback, None),
            '_epics.clear_callback':(_clear_callback, None),
            '_epics.caget':(_caget, None),
            '_epics.caput':(_caput, None),
            '_epics.pend_io':(_pend_io, None),
            '_epics.pend_event':(_pend_event, None) }
