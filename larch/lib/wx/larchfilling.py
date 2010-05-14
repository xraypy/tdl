"""
Larch Filling:  stolen and hacked from PyCrust's Filling module

Filling is the gui tree control through which a user can navigate
the local namespace or any object."""

__author__ = "Patrick K. O'Brien <pobrien@orbtech.com>"
__cvsid__ = "$Id: filling.py 37633 2006-02-18 21:40:57Z RD $"
__revision__ = "$Revision: 37633 $"[11:-2]

import wx
import sys
import types

from wx.py import dispatcher
from wx.py import editwindow

import inspect

from wx.py import introspect
from larch.symboltable import SymbolTable, Group
from larch.closure import Closure

VERSION = '0.9.5(Larch)'

COMMONTYPES = [getattr(types, t) for t in dir(types) \
               if not t.startswith('_') \
               and t not in ('ClassType', 'InstanceType', 'ModuleType')]

DOCTYPES = ('BuiltinFunctionType', 'BuiltinMethodType', 'ClassType',
            'FunctionType', 'GeneratorType', 'InstanceType',
            'LambdaType', 'MethodType', 'ModuleType',
            'UnboundMethodType', 'method-wrapper')

SIMPLETYPES = [getattr(types, t) for t in dir(types) \
               if not t.startswith('_') and t not in DOCTYPES]

del t

try:
    COMMONTYPES.append(type(''.__repr__))  # Method-wrapper in version 2.2.x.
except AttributeError:
    pass

class FillingTree(wx.TreeCtrl):
    """FillingTree based on TreeCtrl."""
    
    name = 'Larch Filling Tree'
    revision = __revision__

    def __init__(self, parent, id=-1, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=wx.TR_DEFAULT_STYLE,
                 rootObject=None, rootLabel=None, rootIsNamespace=False,
                 static=False):
        """Create FillingTree instance."""
        wx.TreeCtrl.__init__(self, parent, id, pos, size, style)
        self.rootIsNamespace = rootIsNamespace
        self.rootLabel = rootLabel
        self.static = static
        self.item = None
        self.root = None
        self.setRootObject(rootObject)

    def setRootObject(self, rootObject=None):
        self.rootObject = rootObject
        if self.rootObject is None:
            return
        if not self.rootLabel:
            self.rootLabel = 'Larch Data'
        rootData = wx.TreeItemData(rootObject)
        self.item = self.root = self.AddRoot(self.rootLabel, -1, -1,  rootData)

        self.SetItemHasChildren(self.root,  self.objHasChildren(self.rootObject))
        self.Bind(wx.EVT_TREE_ITEM_EXPANDING, self.OnItemExpanding, id=self.GetId())
        self.Bind(wx.EVT_TREE_ITEM_COLLAPSED, self.OnItemCollapsed, id=self.GetId())
        self.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnSelChanged, id=self.GetId())
        self.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.OnItemActivated, id=self.GetId())
        self.Bind(wx.EVT_RIGHT_DOWN, self.OnSelChanged, id=self.GetId() )        
        if not self.static:
            dispatcher.connect(receiver=self.push, signal='Interpreter.push')

    def push(self, command, more):
        """Receiver for Interpreter.push signal."""
        self.display()

    def OnItemExpanding(self, event):
        """Add children to the item."""
        busy = wx.BusyCursor()
        item = event.GetItem()
        if self.IsExpanded(item):
            return
        self.addChildren(item)
#        self.SelectItem(item)

    def OnItemCollapsed(self, event):
        """Remove all children from the item."""
        busy = wx.BusyCursor()
        item = event.GetItem()

    def OnSelChanged(self, event):
        """Display information about the item."""
        busy = wx.BusyCursor()
        if hasattr(event, 'GetItem'):
            self.item = event.GetItem()
        self.display()

    def OnItemActivated(self, event):
        """Launch a DirFrame."""
        item = event.GetItem()
        text = self.getFullName(item)
        obj = self.GetPyData(item)
        frame = FillingFrame(parent=self, size=(500, 500),
                             rootObject=obj,
                             rootLabel=text, rootIsNamespace=False)
        frame.Show()

    def objHasChildren(self, obj):
        """Return true if object has children."""
        if self.objGetChildren(obj):
            return True
        else:
            return False

    def objGetChildren(self, obj):
        """Return dictionary with attributes or contents of object."""
        busy = wx.BusyCursor()
        otype = type(obj)
        d = {}
        self.ntop = 0
        if isinstance(obj, SymbolTable)   or isinstance(obj, Group):
            return obj._publicmembers()
        if isinstance(obj, (int, float, bool)) or (isinstance(obj, dict) and hasattr(obj, 'keys')):
            return obj
        elif isinstance(obj, (list, tuple)):
            for n in range(len(obj)):
                key = '[' + str(n) + ']'
                d[key] = obj[n]
            
        elif not hasattr(obj, '__call__'):
            for key in introspect.getAttributeNames(obj):
                # Believe it or not, some attributes can disappear,
                # such as the exc_traceback attribute of the sys
                # module. So this is nested in a try block.
                try:
                    if not (key.startswith('__') and key.endswith('__')):
                        d[key] = getattr(obj, key)
                except:
                    pass
        return d

    def addChildren(self, item):
        self.DeleteChildren(item)
        obj = self.GetPyData(item)
        children = self.objGetChildren(obj)
        if not children:
            return
        try:
            keys = children.keys()
        except:
            return
        keys.sort(lambda x, y: cmp(str(x).lower(), str(y).lower()))
        for key in keys:
            itemtext = str(key)
            # Show string dictionary items with single quotes, except
            # for the first level of items, if they represent a
            # namespace.
            if (type(obj) is types.DictType \
                and type(key) is types.StringType \
                and (item != self.root \
                     or (item == self.root and not self.rootIsNamespace))):
                itemtext = repr(key)
            child = children[key]
            data = wx.TreeItemData(child)
            branch = self.AppendItem(parent=item, text=itemtext, data=data)
            self.SetItemHasChildren(branch, self.objHasChildren(child))

    def display(self):
        item = self.item
        if not item:
            return
        obj = self.GetPyData(item)
        if isinstance(obj, Closure):
            obj = obj.func
            
        if self.IsExpanded(item):
            self.addChildren(item)
        self.setText('')

        if wx.Platform == '__WXMSW__':
            if obj is None: # Windows bug fix.
                return
        self.SetItemHasChildren(item, self.objHasChildren(obj))
        otype = type(obj)
        text = ''
        text += self.getFullName(item)
        text += '\n\nType: ' + str(otype)

        try:
            value = str(obj)
        except:
            value = ''
        if otype is types.StringType or otype is types.UnicodeType:
            value = repr(obj)
        text += '\n\nValue: ' + value
        need_doc = True
        try:
            text += '\n\nDocstring:\n\n"""' + \
                    obj.__doc__.strip() + '"""'
            need_doc = False
        except:
            try:
                text += '\n\nDocstring:\n\n"""' + \
                        obj.__doc__().strip() + '"""'
                need_doc = False
            except:
                pass
            
        if need_doc:
            if otype is types.InstanceType:
                try:
                    text += '\n\nClass Definition:\n\n' + \
                            inspect.getsource(obj.__class__)
                except:
                    pass
            else:
                try:
                    text += '\n\nSource Code:\n\n' + \
                            inspect.getsource(obj)
                except:
                    pass
        self.setText(text)

    def getFullName(self, item, partial=''):
        """Return a syntactically proper name for item."""
        name = self.GetItemText(item)
        parent = None
        obj = None
        if item != self.root:
            parent = self.GetItemParent(item)
            obj = self.GetPyData(parent)
        # Apply dictionary syntax to dictionary items, except the root
        # and first level children of a namepace.
        if (type(obj) is types.DictType \
            or str(type(obj))[17:23] == 'BTrees' \
            and hasattr(obj, 'keys')) \
        and ((item != self.root and parent != self.root) \
            or (parent == self.root and not self.rootIsNamespace)):
            name = '[' + name + ']'
        # Apply dot syntax to multipart names.
        if partial:
            if partial[0] == '[':
                name += partial
            else:
                name += '.' + partial
        # Repeat for everything but the root item
        # and first level children of a namespace.
        if (item != self.root and parent != self.root) \
        or (parent == self.root and not self.rootIsNamespace):
            name = self.getFullName(parent, partial=name)
        return name

    def setText(self, text):
        """Display information about the current selection."""

        # This method will likely be replaced by the enclosing app to
        # do something more interesting, like write to a text control.
        print text

    def setStatusText(self, text):
        """Display status information."""

        # This method will likely be replaced by the enclosing app to
        # do something more interesting, like write to a status bar.
        print text


class FillingText(editwindow.EditWindow):
    """FillingText based on StyledTextCtrl."""

    name = 'Filling Text'
    revision = __revision__

    def __init__(self, parent, id=-1, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=wx.CLIP_CHILDREN,
                 static=False):
        """Create FillingText instance."""
        editwindow.EditWindow.__init__(self, parent, id, pos, size, style)
        # Configure various defaults and user preferences.
        self.SetReadOnly(True)
        self.SetWrapMode(True)
        self.SetMarginWidth(1, 0)
        if not static:
            dispatcher.connect(receiver=self.push, signal='Interpreter.push')

    def push(self, command, more):
        """Receiver for Interpreter.push signal."""
        self.Refresh()

    def SetText(self, *args, **kwds):
        self.SetReadOnly(False)
        editwindow.EditWindow.SetText(self, *args, **kwds)
        self.SetReadOnly(True)


class Filling(wx.SplitterWindow):
    """Filling based on wxSplitterWindow."""

    name = 'Filling'
    revision = __revision__

    def __init__(self, parent, id=-1, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=wx.SP_3D|wx.SP_LIVE_UPDATE,
                 name='Filling Window', rootObject=None,
                 rootLabel=None, rootIsNamespace=False, static=False):
        """Create a Filling instance."""
        wx.SplitterWindow.__init__(self, parent, id, pos, size, style, name)

        self.tree = FillingTree(parent=self, rootObject=rootObject,
                                rootLabel=rootLabel,
                                rootIsNamespace=rootIsNamespace,
                                static=static)
        self.text = FillingText(parent=self, static=static)
        
        wx.FutureCall(1, self.SplitVertically, self.tree, self.text, 200)
        
        self.SetMinimumPaneSize(1)

        # Override the filling so that descriptions go to FillingText.
        self.tree.setText = self.text.SetText

        # Display the root item.
        if self.tree.root is not None:
            self.tree.SelectItem(self.tree.root)
            self.tree.display()

        self.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGED, self.OnChanged)

    def SetRootObject(self, rootObject=None):
        self.tree.setRootObject(rootObject)
        if self.tree.root is not None:        
            self.tree.SelectItem(self.tree.root)
            self.tree.display()
        
    def OnChanged(self, event):
        #this is important: do not evaluate this event=> otherwise, splitterwindow behaves strange
        #event.Skip()
        pass


    def LoadSettings(self, config):
        pos = config.ReadInt('Sash/FillingPos', 200)
        wx.FutureCall(250, self.SetSashPosition, pos)
        zoom = config.ReadInt('View/Zoom/Filling', -99)
        if zoom != -99:
            self.text.SetZoom(zoom)

    def SaveSettings(self, config):
        config.WriteInt('Sash/FillingPos', self.GetSashPosition())
        config.WriteInt('View/Zoom/Filling', self.text.GetZoom())



class FillingFrame(wx.Frame):
    """Frame containing the namespace tree component."""

    name = 'Filling Frame'
    revision = __revision__
    def __init__(self, parent=None, id=-1, title='Larch Data Tree',
                 pos=wx.DefaultPosition, size=(600, 400),
                 style=wx.DEFAULT_FRAME_STYLE, rootObject=None,
                 rootLabel=None, rootIsNamespace=False, static=False):
        """Create FillingFrame instance."""
        wx.Frame.__init__(self, parent, id, title, pos, size, style)
        intro = 'Larch Data Tree (borrowed from PyFilling)'
        self.CreateStatusBar()
        self.SetStatusText(intro)
#         import images
#         self.SetIcon(images.getPyIcon())
        self.filling = Filling(parent=self, rootObject=rootObject,
                               rootLabel=rootLabel,
                               rootIsNamespace=rootIsNamespace,
                               static=static)
        # Override so that status messages go to the status bar.
        self.filling.tree.setStatusText = self.SetStatusText


class App(wx.App):
    """PyFilling standalone application."""
    def OnInit(self):
        wx.InitAllImageHandlers()
        self.fillingFrame = FillingFrame()
        self.fillingFrame.Show(True)
        self.SetTopWindow(self.fillingFrame)
        return True
