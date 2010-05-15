import wx
import time
import os

MODNAME = '_shell'

def ensuremod(larch):
    if larch is not None:
        symtable = larch.symtable
        if not symtable.has_group(MODNAME):
            symtable.newgroup(MODNAME)
    
def _WXcd(parent=None, larch=None, **kws):
    """Directory Browser to Change Directory"""
    dlg = wx.DirDialog(self, message='Choose Directory',
                       defaultDir = os.getcwd())
    if dlg.ShowModal() == wx.ID_OK:
        path = dlg.GetPath()
        os.chdir(path)
    return os.getcwd()

def _WXsave(parent=None, larch=None,
            message = 'Save File As', 
            fname=None, choices=None, **k):
    symtable = ensuremod(larch)
    if fname is None:
        try:
            fname = symtable.get_symbol("%s.deffile" % MODNAME)
        except:
            fname = ''
    if choices  is None:
        try:
            choices = symtable.get_symbol("%s.fchoices" % MODNAME)
        except:
            choices = []
    if '*' not in choices:
        choices.append('*')
    
    dlg = wx.FileDialog(self, message=message, 
                        defaultDir = os.getcwd(),
                        defaultFile= fname,
                        wildcard= choices,
                        style=wx.SAVE|wx.CHANGE_DIR)

    if dlg.ShowModal() == wx.ID_OK:
        path = dlg.GetPath()
        return path
    
def registerPlugin():
    return ('_shell', True,
            {'gcd': _WXcd,
             'savefile': _WXsave})


        
