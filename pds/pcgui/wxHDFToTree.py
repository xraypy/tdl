import wx
import h5py
import time

def dictToTree(thisDict, thisTree, thisRoot):
        for key, value in thisDict.items():
            if type(value) == dict:
                parentItem = thisTree.AppendItem(thisRoot, key)
                thisTree.SetPyData(parentItem, None)
                dictToTree(value, thisTree, parentItem)
            else:
                childItem = thisTree.AppendItem(thisRoot, str(key))
                thisTree.SetPyData(childItem, thisDict[key])
        thisTree.SortChildren(thisRoot)

class hdfToTree(wx.Frame):
    def __init__(self, *args, **kwargs):
        #print time.ctime()
        wx.Frame.__init__(self, args[0], -1, title='Test', size=(200, 400))
        self.fullWindow = wx.Panel(self)
        testFile = \
            h5py.File('C:\\Users\\biwer\\Desktop\\HDFFiles\\ProjectALL.h5', 'r')
        #print len(testFile.items())
        self.treeSizer = wx.BoxSizer(wx.VERTICAL)
        self.scanTree = myTreeCtrl(self.fullWindow)
        self.testRoot = self.scanTree.AddRoot('testing')
        testDict = {}
        for item in testFile.items():
            item0 = item[0]
            item1 = item[1]
            itemString = item1.attrs.get('name')
            (specName, scanNum, pointNum, epoch) = itemString.split(':')
            scanNum = 'Scan ' + scanNum[1:]
            pointNum = pointNum.split('/')[0]
            pointNum = 'Point ' + pointNum[1:]
            if specName not in testDict:
                testDict[specName] = {scanNum: {pointNum: item0}}
            elif scanNum not in testDict[specName]:
                testDict[specName][scanNum] = {pointNum: item0}
            elif pointNum not in testDict[specName][scanNum]:
                testDict[specName][scanNum][pointNum] = item0
            else:
                print 'Error'
                break
        dictToTree(testDict, self.scanTree, self.testRoot)
        self.Bind(wx.EVT_TREE_SEL_CHANGED, self.newSelected, self.scanTree)
        self.treeSizer.Add(self.scanTree, proportion=1, flag=wx.EXPAND | wx.ALL)
        self.fullWindow.SetSizerAndFit(self.treeSizer)
        self.Show()
        #print time.ctime()

    def newSelected(self, event):
        ofMe = event.GetItem()
        print self.scanTree.GetItemPyData(ofMe)

class myTreeCtrl(wx.TreeCtrl):
    def __init__(self, *args, **kwargs):
        wx.TreeCtrl.__init__(self, *args, **kwargs)
        
    def OnCompareItems(self, item1, item2):
        item1Str = self.GetItemText(item1)
        item2Str = self.GetItemText(item2)
        if item1Str.startswith('Scan') or item1Str.startswith('Point'):
            item1Num = item1Str.split(' ')[1]
            item2Num = item2Str.split(' ')[1]
            return cmp(int(item1Num), int(item2Num))
        else:
            return cmp(item1Str, item2Str)


if __name__ == '__main__':
    app = wx.App()
    testing = hdfToTree(None)
    app.MainLoop()
