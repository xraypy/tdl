import wx
import time

class hdfToTree():
    def __init__(self):
        self.reverseLookup = {}
        
    def dictToTree(self, thisDict, thisTree, thisRoot):
        """Turns a nested dictionary into a tree under thisRoot."""
        for key, value in thisDict.items():
            if type(value) == dict:
                parentItem = thisTree.AppendItem(thisRoot, key)
                thisTree.SetPyData(parentItem, None)
                self.dictToTree(value, thisTree, parentItem)
            else:
                childItem = thisTree.AppendItem(thisRoot, str(key))
                thisTree.SetPyData(childItem, thisDict[key])
        thisTree.SortChildren(thisRoot)
    
    def populateTree(self, thisTree, thisObject):
        """Given a tree and an hdf_data object, build a 
            dictionary out of the object and use it to
            populate the tree using dictToTree."""
        rootName = thisObject.fname.split('/')[-1].split('\\')[-1]
        treeRoot = thisTree.AddRoot(rootName)
        hdfDict = {}
        for item in thisObject.all_items:
            item0 = item[0]
            item1 = item[1]
            itemString = item1.attrs.get('name')
            (specName, scanNum, pointNum, epoch) = itemString.split(':')
            scanNum = 'Scan ' + scanNum[1:]
            pointNum = pointNum.split('/')[0]
            pointNum = 'Point ' + pointNum[1:]
            if specName not in hdfDict:
                hdfDict[specName] = {scanNum: {pointNum: item0}}
            elif scanNum not in hdfDict[specName]:
                hdfDict[specName][scanNum] = {pointNum: item0}
            elif pointNum not in hdfDict[specName][scanNum]:
                hdfDict[specName][scanNum][pointNum] = item0
            else:
                print 'Error: Duplicate point'
                print 'Specfile: ' + specName
                print 'Scan number: ' + scanNum
                print 'Point number: ' + pointNum
                continue
        self.dictToTree(hdfDict, thisTree, treeRoot)
        self.populateReverse(thisTree, treeRoot)
    
    def populateReverse(self, thisTree, thisRoot):
        item, cookie = thisTree.GetFirstChild(thisRoot)
        while item:
            iterData = thisTree.GetItemPyData(item)
            if iterData is not None:
                self.reverseLookup[iterData] = item
            else:
                self.populateReverse(thisTree, item)
            item, cookie = thisTree.GetNextChild(thisRoot, cookie)
    
    def deleteItem(self, thisTree, thisObject, thisItem):
        """Given a tree, an hdf_data object, and an item,
            delete the item and its children from both
            the tree and the object."""
        thisData = thisTree.GetItemPyData(thisItem)
        if thisData is not None:
            thisObject.delete(thisData)
        elif thisTree.ItemHasChildren(thisItem):
            item, cookie = thisTree.GetFirstChild(thisItem)
            while item:
                item2, cookie2 = thisTree.GetNextChild(thisItem, cookie)
                self.deleteItem(thisTree, thisObject, item)
                item, cookie = item2, cookie2
        thisTree.Delete(thisItem)
    
    def getRelevantChildren(self, thisTree, thisObject, thisItem):
        """Given a tree, an hdf_data object, and an item,
            return a list of all the item's children that
            are associated with a point in the object.
            This list includes the point itself to allow
            for recursive function calls."""
        thisData = thisTree.GetItemPyData(thisItem)
        if thisData is not None:
            return [thisData]
        elif thisTree.ItemHasChildren(thisItem):
            returnList = []
            item, cookie = thisTree.GetFirstChild(thisItem)
            while item:
                returnList.extend(self.getRelevantChildren(thisTree,
                                                           thisObject, item))
                item, cookie = thisTree.GetNextChild(thisItem, cookie)
            return returnList
        else:
            return []
    
    def statusString(self, thisTree, thisObject, thisItem):
        """Given a tree, an hdf_data object, and an item,
            return a string describing the item."""
        thisData = thisTree.GetItemPyData(thisItem)
        if thisData is not None:
            thisName = thisObject[thisData]['name']
            thisSpec, thisScan, thisPoint, thisTime = thisName.split(':')
            thisScan = thisScan[1:]
            thisPoint = thisPoint.split('/')[0][1:]
            return 'Scan ' + thisScan + ', Point ' + thisPoint
        else:
            return thisTree.GetItemText(thisItem)

    '''def newSelected(self, event):
        ofMe = event.GetItem()
        print self.scanTree.GetItemPyData(ofMe)'''

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
