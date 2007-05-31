############################################################################
# Tom Trainor (fftpt@uaf.edu)
# This should be run as a child window from wxGUI
#
#
############################################################################

from PythonCard import model, dialog
import wx
import os
import xrf_lookup
from wxTdlUtil import wxTdlUtil
from Util import PrintExceptErr
import numpy as Num

"""
Plot GUI
"""
############################################################################
#
#############################################################################

class wxPlotSelection(model.Background, wxTdlUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.components.PlotParams._autoresize = 0
        self.components.AxisParams._autoresize = 0

        self.startup = True
        
        # set up tdl
        self.shell = None
        self.tdl = None
        self.init_tdl()

        #
        self.components.Line.items = ['-','--','-.',':']
        self.components.Symbol.items = ['.',',','o','^','v','<','>','s','+','x',
                                        'D','d','1','2','3','4','h','H','p','|','_']
        self.components.Color.items = ['b','g','r','c','m','y','k','w']

        # init all the menus
        self.init_tdl_list_items()

    def init_tdl_list_items(self):
        # the below lists are just updated
        # to include tdl variables ie list items
        # (should not change selections)
        self.init_SaveNameItems()
        self.init_XItems()
        self.init_YItems()
        self.init_XerrItems()
        self.init_YerrItems()
        self.init_PlotNumItems()
        self.init_PlotNum_axisItems()

    def updateOnEnter(self,event):
        "Update the list"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_PlotNumItems()
            self.PlotParams_update_from_fields()
        else:
            event.skip()
        return

    def updateAxisOnEnter(self,event):
        "Update the list"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_PlotNum_axisItems()
            self.AxisParams_update_from_fields()
        else:
            event.skip()
        return

    def calc_PlotNumList(self):
        nr = self.components.Nrow.text.strip()
        nc = self.components.Ncol.text.strip()
        sel = self.components.PlotNum.stringSelection
        if len(nr) == 0:
            nr = '1'
            self.components.Nrow.text = nr
        if len(nc) == 0:
            nc = '1'
            self.components.Ncol.text = nc
        nplot = int(nr)*int(nc)
        list = Num.arange(1,nplot+1,dtype='int')
        list = map(str,list)
        return list
    ###########################################################
    #             Menus                                       #
    ###########################################################
    def on_menuFileExit_select(self,event): 
        #self.close()
        self.exit()
        
    ###########################################################
    #       Update Tdl Variables              
    ###########################################################
    def on_UpdateTdl_mouseClick(self, event):
        "use this to force update incase the data list has changed"
        self.init_tdl_list_items()
        self.PlotParams_update_from_fields()
        self.AxisParams_update_from_fields()
        return

    ###########################################################
    #   PlotNum / Ncol / Nrow
    ###########################################################
    #def on_PlotNum_keyDown(self,event):
    #    "Update the list"
    #    self.updateOnEnter(event)
    #    return
    def init_PlotNumItems(self):
        list = self.calc_PlotNumList()
        sel = self.components.PlotNum.stringSelection
        self.components.PlotNum.items = list
        self.components.PlotNum.stringSelection = sel
        return

    def on_PlotNum_select(self,event):
        "Update the list"
        self.PlotParams_update_from_fields()
        return

    def on_Ncol_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        self.init_PlotNum_axisItems()
        return

    def on_Nrow_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        self.init_PlotNum_axisItems()
        return

    ###########################################################
    #   PlotLabel
    ###########################################################
    def on_Label_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    ###########################################################
    #   X/Y
    ###########################################################
    def init_XItems(self):
        " Initialize the menu    "
        X_str = self.components.X.text
        self.components.X.items = self.listAllData()
        self.components.X.text = X_str
        return

    def on_X_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    def init_YItems(self):
        " Initialize the menu    "
        Y_str = self.components.Y.text
        self.components.Y.items = self.listAllData()
        self.components.Y.text = Y_str
        return

    def on_Y_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    ###########################################################
    #   Xerr/Yerr
    ###########################################################
    def init_XerrItems(self):
        " Initialize the menu    "
        Xerr_str = self.components.Xerr.text
        self.components.Xerr.items = self.listAllData()
        self.components.Xerr.text = Xerr_str
        return

    def on_Xerr_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    def init_YerrItems(self):
        " Initialize the menu    "
        Yerr_str = self.components.Yerr.text
        self.components.Yerr.items = self.listAllData()
        self.components.Yerr.text = Yerr_str
        return

    def on_Yerr_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    ###########################################################
    #   X/Yscale
    ###########################################################
    def on_Xscale_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    def on_Yscale_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    ###########################################################
    #   X/Ynorm
    ###########################################################
    #def on_Xnorm_select(self,event):
    #    "Update the list"
    #    self.PlotParams_update_from_fields()
    #    return

    #def on_Ynorm_select(self,event):
    #    "Update the list"
    #    self.PlotParams_update_from_fields()
    #    return

    ###########################################################
    #   X/Yoff
    ###########################################################
    def on_Xoff_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    def on_Yoff_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    ###########################################################
    #   KW
    ###########################################################
    def on_KwArgs_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    ###########################################################
    #   color/line/symbol
    ###########################################################
    def on_Color_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    def on_Line_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return

    def on_Symbol_keyDown(self,event):
        "Update the list"
        self.updateOnEnter(event)
        return
    
    ###########################################################
    #       Plot Parameters                                   
    ###########################################################
    def on_PlotParams_select(self,event):
        #print self.components.PkParams.items
        selected =  self.components.PlotParams.getStringSelection()
        plot_params = selected[0]
        #print plot_params
        self.put_PlotPar_fields(plot_params)
        return

    def on_Delete_mouseClick(self,event):
        selected =  self.components.PlotParams.getStringSelection()
        if selected == None: return
        list = self.components.PlotParams.items
        for sel in selected:
            j = 0
            while j < len(list):
                if sel[0] == list[j][0]:
                    list.pop(j)
                    break
                j = j + 1
        self.components.PlotParams.items = list
        self.components.PlotParams.SetSelection(0)
        return

    def get_PlotParams(self):
        """ return list of all Peak Parameters from GUI """
        return self.components.PlotParams.items

    def PlotParams_update_from_fields(self):
        plot_params = self.get_PlotPar_fields()
        self.PlotParams_update(plot_params)
    
    def PlotParams_update(self,plot_params):        
        "update the pk params list"
        list = self.components.PlotParams.items

        if len(list) == 0:
            list.append(plot_params)
            sel = 0
        else:
            found = False
            for j in range(len(list)):
                if (list[j][0] == plot_params[0]) and (list[j][1] == plot_params[1]):
                    list[j] = plot_params
                    found = True
                    sel = j
                    break
            if found == False:
                list.append(plot_params)
        self.components.PlotParams.items = list
        return

    def get_PlotPar_fields(self):
        "ret a list of all the peak parameter info in the entry fields"
        plot_params = ['']*19
        # gather all the info
        plot_params[0] = self.components.Label.text.strip()        
        plot_params[1] = self.components.PlotNum.stringSelection.strip()
        plot_params[2] = self.components.X.text.strip()
        plot_params[3] = self.components.Y.text.strip()
        plot_params[4] = self.components.Xerr.text.strip()
        plot_params[5] = self.components.Yerr.text.strip()
        plot_params[6] = self.components.Xscale.text.strip()
        plot_params[7] = self.components.Yscale.text.strip()
        plot_params[8] = str(self.components.Xnorm.checked)
        plot_params[9] = str(self.components.Ynorm.checked)
        plot_params[10] = self.components.Xoff.text.strip()
        plot_params[11] = self.components.Yoff.text.strip()
        plot_params[12] = str(self.components.Xcen.checked)
        plot_params[13] = str(self.components.Ycen.checked)
        plot_params[14] = self.components.Color.text.strip()
        plot_params[15] = self.components.Line.text.strip()
        plot_params[16] = self.components.Symbol.text.strip()
        plot_params[17] = self.components.KwArgs.text.strip()
        plot_params[18] = str(self.components.Ignore.checked)
        return plot_params

    def put_PlotPar_fields(self,plot_params):
        "reverse of above"
        if plot_params == None: return
        self.components.Label.text      = plot_params[0]
        self.components.PlotNum.stringSelection    = plot_params[1]
        self.components.X.text          = plot_params[2]
        self.components.Y.text          = plot_params[3] 
        self.components.Xerr.text       = plot_params[4] 
        self.components.Yerr.text       = plot_params[5] 
        self.components.Xscale.text     = plot_params[6] 
        self.components.Yscale.text     = plot_params[7] 
        self.components.Xnorm.checked   = eval(plot_params[8])  
        self.components.Ynorm.checked   = eval(plot_params[9])  
        self.components.Xoff.text       = plot_params[10] 
        self.components.Yoff.text       = plot_params[11] 
        self.components.Xcen.checked    = eval(plot_params[12]) 
        self.components.Ycen.checked    = eval(plot_params[13])
        self.components.Color.text      = plot_params[14] 
        self.components.Line.text       = plot_params[15]
        self.components.Symbol.text     = plot_params[16]
        self.components.KwArgs.text     = plot_params[17]
        self.components.Ignore.checked  = eval(plot_params[18])
        return

    ###########################################################
    # Axis
    ###########################################################

    #def on_PlotNum_axis_keyDown(self,event):
    #    "Update the list"
    #    self.updateAxisOnEnter(event)
    #    return
    def init_PlotNum_axisItems(self):
        list = self.calc_PlotNumList()
        sel = self.components.PlotNum_axis.stringSelection
        self.components.PlotNum_axis.items = list
        self.components.PlotNum_axis.stringSelection = sel
        return

    def on_PlotNum_axis_select(self,event):
        "Update the list"
        self.AxisParams_update_from_fields()
        return

    def on_Xmin_keyDown(self,event):
        "Update the list"
        self.updateAxisOnEnter(event)
        return

    def on_Xmax_keyDown(self,event):
        "Update the list"
        self.updateAxisOnEnter(event)
        return

    def on_Ymin_keyDown(self,event):
        "Update the list"
        self.updateAxisOnEnter(event)
        return

    def on_Ymax_keyDown(self,event):
        "Update the list"
        self.updateAxisOnEnter(event)
        return

    def on_Xtitle_keyDown(self,event):
        "Update the list"
        self.updateAxisOnEnter(event)
        return

    def on_Ytitle_keyDown(self,event):
        "Update the list"
        self.updateAxisOnEnter(event)
        return

    def on_PlotTitle_keyDown(self,event):
        "Update the list"
        self.updateAxisOnEnter(event)
        return

    def on_LegendLoc_keyDown(self,event):
        "Update the list"
        self.updateAxisOnEnter(event)
        return

    ########
    def on_AxisParams_select(self,event):
        #print self.components.PkParams.items
        selected =  self.components.AxisParams.getStringSelection()
        axis_params = selected[0]
        #print plot_params
        self.put_AxisPar_fields(axis_params)
        return

    def on_DeleteAxis_mouseClick(self,event):
        selected =  self.components.AxisParams.getStringSelection()
        if selected == None: return
        list = self.components.AxisParams.items
        for sel in selected:
            j = 0
            while j < len(list):
                if sel[0] == list[j][0]:
                    list.pop(j)
                    break
                j = j + 1
        self.components.AxisParams.items = list
        self.components.AxisParams.SetSelection(0)
        return

    def get_AxisParams(self):
        """ return list of all Peak Parameters from GUI """
        return self.components.AxisParams.items

    def AxisParams_update_from_fields(self):
        axis_params = self.get_AxisPar_fields()
        self.AxisParams_update(axis_params)
    
    def AxisParams_update(self,axis_params):        
        "update the axis params list"
        list = self.components.AxisParams.items

        if len(list) == 0:
            list.append(axis_params)
            sel = 0
        else:
            found = False
            for j in range(len(list)):
                if list[j][0] == axis_params[0]:
                    list[j] = axis_params
                    found = True
                    sel = j
                    break
            if found == False:
                list.append(axis_params)
        self.components.AxisParams.items = list
        return

    def get_AxisPar_fields(self):
        "ret a list of all the peak parameter info in the entry fields"
        axis_params = ['']*13
        # gather all the info
        axis_params[0] = self.components.PlotNum_axis.stringSelection.strip()        
        axis_params[1] = self.components.Xmin.text.strip()
        axis_params[2] = self.components.Xmax.text.strip()
        axis_params[3] = self.components.Ymin.text.strip()
        axis_params[4] = self.components.Ymax.text.strip()
        axis_params[5] = str(self.components.Xlog.checked)
        axis_params[6] = str(self.components.Ylog.checked)
        axis_params[7] = self.components.Xtitle.text.strip()
        axis_params[8] = self.components.Ytitle.text.strip()
        axis_params[9] = self.components.PlotTitle.text.strip()
        axis_params[10] = str(self.components.ShowLegend.checked)
        axis_params[11] = self.components.LegendLoc.text.strip()
        axis_params[12] = str(self.components.AxisIgnore.checked)
        return axis_params

    def put_AxisPar_fields(self,axis_params):
        "reverse of above"
        if axis_params == None: return
        self.components.PlotNum_axis.stringSelection = axis_params[0]
        self.components.Xmin.text         = axis_params[1]
        self.components.Xmax.text         = axis_params[2]
        self.components.Ymin.text         = axis_params[3] 
        self.components.Ymax.text         = axis_params[4] 
        self.components.Xlog.checked      = eval(axis_params[5]) 
        self.components.Ylog.checked      = eval(axis_params[6]) 
        self.components.Xtitle.text       = axis_params[7] 
        self.components.Ytitle.text       = axis_params[8]  
        self.components.PlotTitle.text    = axis_params[9]  
        self.components.ShowLegend.checked  = eval(axis_params[10]) 
        self.components.LegendLoc.text    = axis_params[11]
        self.components.AxisIgnore.checked  = eval(axis_params[12]) 
        return

    ###########################################################
    # Save/Restore
    ###########################################################
    def init_SaveNameItems(self):
        " Initialize the menu    "
        SaveName_str = self.components.SaveName.text
        self.components.SaveName.items = self.listAllData()
        self.components.SaveName.text = SaveName_str
        return

    def on_Save_mouseClick(self,event):
        save_name = self.components.SaveName.text.strip()
        if len(save_name) == 0:
            self.components.PrintCmds.appendText('No Save Name\n')
            return
        plot_params = self.get_PlotParams()
        axis_params = self.get_AxisParams()
        nr = self.components.Nrow.text.strip()
        nc = self.components.Ncol.text.strip()
        if len(nr) == 0: nr = '1'
        if len(nc) == 0: nc = '1'
        data = {'plot_params':plot_params,'axis_params':axis_params,'nr':nr,'nc':nc}
        self.setValue(save_name,data)
        return

    def on_Load_mouseClick(self,event):
        save_name = self.components.SaveName.text.strip()
        if len(save_name) == 0:
            self.components.PrintCmds.appendText('No Save Name\n')
            return
        data = self.getValue(save_name)
        try:
            plot_params = data['plot_params']
            self.components.PlotParams.items = plot_params
            axis_params = data['axis_params']
            self.components.AxisParams.items = axis_params
            nr = data['nr']
            self.components.Nrow.text = nr
            nc = data['nc']
            self.components.Ncol.text = nc
        except:
            self.components.PrintCmds.appendText('Failed to load plot parameters\n')
        #
        self.components.PlotParams.SetSelection(0)
        self.components.AxisParams.SetSelection(0)
        return

    ###########################################################
    # Plotting
    ###########################################################
    # note, Plot function should pack up data similiar to the save
    # function, and we can then move plot_cmd to the Pylab module
    # and call from here
    def on_Plot_mouseClick(self,event):
        "make a plot"
        # update data and plot params
        self.PlotParams_update_from_fields()
        self.AxisParams_update_from_fields()
        self.plot_cmd()
        return

    def plot_cmd(self):
        self.execLine("pylab.clf()")
        #self.execLine("pylab.hold(True)")
        nr = self.components.Nrow.text.strip()
        nc = self.components.Ncol.text.strip()
        if len(nr) == 0: nr = '1'
        if len(nc) == 0: nc = '1'

        plot_params = self.components.PlotParams.items
        counter = 0
        for params in plot_params:
            Label =  params[0]
            plotnum  = params[1]
            X     =  params[2]
            Y     =  params[3] 
            Xerr  =  params[4] 
            Yerr  =  params[5] 
            Xscale = params[6] 
            Yscale = params[7] 
            Xnorm  = params[8]  
            Ynorm  = params[9]  
            Xoff   = params[10] 
            Yoff   = params[11] 
            Xcen   = params[12] 
            Ycen   = params[13]
            color  = params[14]
            line   = params[15]
            symbol = params[16]
            kw      = params[17]
            ignore  = eval(params[18])
            if ignore != True:
                #build and exec plot str
                #plotter(x,y=None,fmt='k-',xerr=None,yerr=None,xscale=1,yscale=1,xnorm=False,ynorm=False,
                #        xoff=0,yoff=0,xcen=False,ycen=False,xlog=False,ylog=False,nr=1,nc=1,np=1,hold=True,**kw)
                if len(X) == 0:
                    pstr = "pylab.plotter(%s," % (Y)
                elif len(Y) == 0:
                    pstr = "pylab.plotter(%s," % (X)
                else:
                    pstr = "pylab.plotter(%s,%s," % (X,Y)
                pstr = pstr + 'fmt=' + "'" + color+symbol+line + "'" + ','
                if len(Xerr) > 0:
                    pstr = pstr + 'xerr=' + Xerr + ','
                if len(Yerr) > 0:
                    pstr = pstr + 'yerr=' + Yerr + ','
                if len(Xscale) > 0:
                    pstr = pstr + 'xscale=' + Xscale + ','
                if len(Yscale) > 0:
                    pstr = pstr + 'yscale=' + Yscale + ','
                pstr = pstr + 'xnorm=' + str(Xnorm) + ','
                pstr = pstr + 'ynorm=' + str(Ynorm) + ','
                if len(Xoff) > 0:
                    pstr = pstr + 'xoff=' + Xoff + ','
                if len(Yoff) > 0:
                    pstr = pstr + 'yoff=' + Yoff + ','
                pstr = pstr + 'xcen=' + str(Xcen) + ','
                pstr = pstr + 'ycen=' + str(Ycen) + ','
                #pstr = pstr + 'xlog=' + str(Xlog) + ','
                #pstr = pstr + 'ylog=' + str(Ylog) + ','
                pstr = pstr + 'nr=' + nr + ','
                pstr = pstr + 'nc=' + nc + ','
                if len(plotnum) == 0: plotnum = '1'
                pstr = pstr + 'np=' + plotnum + ','
                pstr = pstr + 'hold=True,'
                if len(Label) == 0: Label = 'Plot %d' % counter
                pstr = pstr + 'label="%s"' % (Label)
                if len(kw) > 0:
                    pstr = pstr + ',' + kw + ')'
                else:
                    pstr = pstr + ')'

                # show print cmd
                self.components.PrintCmds.appendText(pstr + '\n')
                # execute
                try:
                    self.execLine(pstr)
                except:
                    PrintExceptErr('Plot command failed')
                counter = counter + 1
        self.components.PrintCmds.appendText('----\n')

        # now do the axis
        axis_params = self.components.AxisParams.items
        for axis in axis_params:
            PlotNum      = axis[0]
            Xmin         = axis[1]
            Xmax         = axis[2]
            Ymin         = axis[3] 
            Ymax         = axis[4] 
            Xlog         = eval(axis[5]) 
            Ylog         = eval(axis[6]) 
            Xtitle       = axis[7] 
            Ytitle       = axis[8]  
            PlotTitle    = axis[9]  
            ShowLegend   = eval(axis[10])
            LegendLoc    = axis[11]
            Ignore       = eval(axis[12])
            if Ignore != True:
                #select plot
                if len(PlotNum) == 0: PlotNum = '1'
                pstr = "pylab.subplot(%s,%s,%s)" % (nr,nc,PlotNum)
                self.execLine(pstr)
                #axis 
                if len(Xmin) or len(Xmin) or len(Xmin) or len(Xmin):
                    #use these as defaults
                    self.execLine("pylab.axis_lims = pylab.axis()")
                    def_axis = self.getValue("pylab.axis_lims")
                    if len(Xmin) == 0: Xmin = def_axis[0]
                    if len(Xmax) == 0: Xmax = def_axis[1]
                    if len(Ymin) == 0: Ymin = def_axis[2]
                    if len(Ymax) == 0: Ymax = def_axis[3]
                    pstr = "pylab.axis([%s,%s,%s,%s])" % (Xmin,Xmax,Ymin,Ymax)
                    self.execLine(pstr)
                # logs
                if Xlog:
                    self.execLine("pylab.semilogx()")
                if Ylog:
                    self.execLine("pylab.semilogy()")
                # titles
                if len(Xtitle) > 0:
                    pstr = 'pylab.xlabel("%s")' % Xtitle
                    self.execLine(pstr)
                if len(Ytitle) > 0:
                    pstr = 'pylab.ylabel("%s")' % Ytitle
                    self.execLine(pstr)
                if len(PlotTitle) > 0:
                    pstr = 'pylab.title("%s")' % PlotTitle
                    self.execLine(pstr)
                # legend
                if ShowLegend:
                    pstr = "pylab.legend("
                    if len(LegendLoc) > 0:
                        pstr = pstr + 'loc=%s' % (LegendLoc)
                    pstr = pstr + ")"
                    self.execLine(pstr)
        return

##################################################
if __name__ == '__main__':
    app = model.Application(wxPlotSelection)
    app.MainLoop()
