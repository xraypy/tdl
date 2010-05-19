"""
Support for Excel files

Notes:
------
This module requires win32com module

Note see pyExcelerator at source forge for alternative


"""
#################################################################################

import numpy as num
import types
import os

#################################################################################
def excel_data(fname,rows,cols,sheet=1, ret_fmt=None,close=True):
    """
    read an excel file
    
    data = excel_data(fname,rows,cols,sheet=1)
    rows = [start_row,end_row],
    cols = [start_col,end_col]
    ret == None: default return (mixed str and float)
    ret == 'f': ret as all floats, strings and Nones are zeros
    ret == 's': ret as all as strings
    THIS IS PRETTY SLOW (at least if excel is not already running?)
    """

    fname = os.path.abspath(fname)
    xls = Excel.Excel(fname)
    if xls == None: return
    
    # if single val make tuple
    if (type(rows) == types.IntType) or (type(rows) == types.FloatType):
        rows = (rows,rows)
    if (type(cols) == types.IntType) or (type(cols) == types.FloatType):
        rows = (rows,rows)
    # make sure we're a tuple
    rows = tuple(rows)
    cols = tuple(cols)
    if (len(rows) != 2) or len(cols)!= 2:
        return
    data = xls.get_range(sheet, rows[0], cols[0], rows[1], cols[1])

    if close:
        xls.close()

    # format ret data
    if ret_fmt == None:
        ret_data = []
        for j in range(len(data)):
            ret_data.append(list(data[j]))
        return ret_data
    elif ret_fmt == 'f':
        try:
            data = num.array(data,dtype=float)
            return data
        except:
            print 'doing float loop'
            ret_data = []
            for j in range(len(data)):
                rdata = []
                for k in range(len(data[j])):
                    try:
                        rdata.append(float(data[j][k]))
                    except:
                        rdata.append(0.0)
                ret_data.append(rdata)
            return num.array(ret_data,dtype=float)
    elif ret_fmt == 's':
        #try:
        #    data = num.array(data,dtype=str)
        #    return data
        #except:
        #print 'doing string loop'
        ret_data = []
        for j in range(len(data)):
            rdata = []
            for k in range(len(data[j])):
                try:
                    rdata.append(str(data[j][k]))
                except:
                    rdata.append('')
            ret_data.append(rdata)
        #return num.array(ret_data,dtype=str)
        return ret_data

    return None

#################################################################################

import win32com.client
import win32com.client.dynamic
import os

class Excel:
    """
    Excel file class
    
    This is Derived from:
    'Python Programming on Win32' by Mark Hammond and Andy Robinson
    """
    def __init__(self, filename=None):
        self.xlApp = win32com.client.dynamic.Dispatch('Excel.Application')
        if filename:
            filename = os.path.abspath(filename)
            self.filename = filename
            self.xlBook = self.xlApp.Workbooks.Open(filename)
        else:
            self.xlBook = self.xlApp.Workbooks.Add()
            self.filename = ''

    def save(self, newfilename=None):
        if newfilename:
            self.filename = newfilename
            self.xlBook.SaveAs(newfilename)
        else:
            self.xlBook.Save()

    def close(self):
        self.xlBook.Close(SaveChanges=0)
        del self.xlApp

    def show(self):
        self.xlApp.Visible = 1

    def hide(self):
        self.xlApp.Visible = 0

    def get_cell(self, sheet, row, col):
        "get value of one cell"
        sht = self.xlBook.Worksheets(sheet)
        return sht.Cells(row, col).Value

    def set_cell(self, sheet, row, col, value):
        "set value of one cell"
        sht = self.xlBook.Worksheets(sheet)
        sht.Cells(row, col).Value = value

    def get_range(self, sheet, row1, col1, row2, col2):
        "return a 2d array (i.e. tuple of tuples)"
        sht = self.xlBook.Worksheets(sheet)
        return sht.Range(sht.Cells(row1, col1), sht.Cells(row2,col2)).Value

    def set_range(self, sheet, leftCol, topRow, data):
        bottomRow = topRow + len(data) - 1
        rightCol = leftCol + len(data[0]) - 1
        sht = self.xlBook.Worksheets(sheet)
        sht.Range(sht.Cells(topRow, leftCol), sht.Cells(bottomRow,rightCol)).Value = data

#################################################################################
if __name__ == "__main__":
    import os
    fname = os.path.join(os.getcwd(),'test.xls')
    xls = Excel(fname)
    sr = 1  #start row
    er = 9  #end row
    sc = 1  #start col
    ec = 2  #end col
    sheet = 1
    data = xls.get_range(sheet,sr,sc,er,ec)
    for j in range(len(data)): print data[j]
    
