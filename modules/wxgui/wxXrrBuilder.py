########################################################################
"""
Tom Trainor (fftpt@uaf.edu)
This should be run as a child window from wxGUI

Modifications:
--------------

"""
########################################################################

from PythonCard import model, dialog
import wx
import os
import types
import copy
import time
import numpy as Num

from   wxUtil import wxUtil
import compound
import interface_model

#########################################################################

class wxXrrBuilder(model.Background, wxUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.components.LayerList._autoresize = 0
        self.startup  = True
        self.dir      = '.'

        # set up shell
        self.shell = None
        self.init_shell()

        # Make sure interface model is loaded
        self.exec_line("import interface_model")

        # init all the menus
        self.init_ShellItems()
        self.init_GUI()
        
    ######################################################

    ###########################################################
    #   Init/Save
    ###########################################################

    ######################################################

    def init_ShellItems(self,):
        self.init_Grp()
        self.init_Model()
        self.init_Data()

    def on_UpdateShell_mouseClick(self,event):
        self.init_ShellItems()

    def on_NewModel_mouseClick(self,event):
        self.components.Grp.text = ''
        self.components.Model.text = ''
        self.components.SaveName.text = ''
        self.components.SlabDelta.text = '10.0'
        self.components.Theta.text = 'Num.arange(0.01,1.01,0.01)'
        self.init_ShellItems()
        self.init_GUI()

    def on_SaveModel_mouseClick(self,event):
        self.UpdateModelFromGui()
    
    def init_GUI(self):
        if self.startup:
            self.startup = False
        
        self.init_LayerList()
        self.init_LayerParams()
        self.init_Params()

        check = self.CheckModelVar()
        if check == False:
            return
        else:
            self.UpdateGuiFromModel()

    ######################################################

    ###########################################################
    #   Menu
    ###########################################################

    ######################################################
    def on_menuFileExit_select(self,event): 
        self.close()

    def on_menuHelpParams_select(self,event): 
        import wxXrrBuilderHelp
        wxXrrBuilderHelp = mod_import(wxXrrBuilderHelp)
        dir       = os.path.dirname(wxXrrBuilderHelp.__file__)
        filename  = os.path.join(dir,'wxXrrBuilderHelp.rsrc.py')
        wxXrrBuilderHelp = wxXrrBuilderHelp.wxXrrBuilderHelp
        self.wxXrrBuilderHelp = model.childWindow(self,wxXrrBuilderHelp,
                                                  filename=filename)
        self.wxXrrBuilderHelp.position = (200, 5)
        self.wxXrrBuilderHelp.visible = True

    ######################################################

    ###########################################################
    #   Update Model/GUI
    ###########################################################

    ######################################################
    def UpdateGuiFromModel(self,):
        """
        update gui from a model
        """
        
        model = self.get_Model()
        if model == None: return

        # layers
        ldata = []        
        for j in range(len(model.layer)):
            ll = []
            m  = model.layer[j]
            ll.append(str(j))
            ll.append(m.name)
            ll.append(str(m.thickness))
            ll.append(str(m.density))
            ll.append(str(m.roughness))
            ldata.append(ll)
        ldata = self.sort_LayerList(ldata)
        self.components.LayerList.items = ldata

        # slab delta
        if model.slab:
            self.components.SlabDelta.text = str(model.slab.delta)

        # Data
        self.init_Data()
        if self.components.Theta.text =='':
            tname = self.get_ModelName()
            tname = tname + '.theta'
            self.components.Theta.text = tname

        # parameters
        params = model.get_params()
        if len(params) > 0:
            self.components.Energy.text    = str(params['energy']) 
            self.components.ConvWidth.text = str(params['wconv'])
            self.components.SampLen.text   = str(params['slen'])
            self.components.BeamVert.text  = str(params['bvert'])
            self.components.BeamHorz.text  = str(params['bhorz'])
            self.components.AreaFlag.stringSelection = str(params['aflag'])
            self.components.RefScale.text  = str(params['rscale'])
            #self.components.FyIndex.text   = str(params['fyidx'])
            #
            if params.has_key('fyel'):
                self.components.FyEl.text = str(params['fyel'])
            else:
                idx = int(params['fyidx'])
                if idx == -1:
                    self.components.FyEl.text   = str(0.0)
                else:
                    elz = model.slab.list_elements()
                    if idx in range(len(elz)):
                        self.components.FyEl.text   = str(elz[idx])
                    else:
                        self.components.FyEl.text   = str(0.0)
            #
            self.components.FyEnergy.text  = str(params['fyenergy'])
            self.components.DetAngle.text  = str(params['adet'])
            self.components.ThetaNorm.text = str(params['tnorm'])
            self.components.RoughFlag.stringSelection = str(params['rflag'])
            self.components.DelZ.text      = str(params['delz'])
            self.components.PDepth.text    = str(params['pdepth'])

    ######################################################
    def UpdateModelFromGui(self):
        """
        create a new model instance
        """
        mname = self.components.SaveName.text.strip()
        if len(mname) == 0:
            print "Model save name required"
            return
        #
        params = self.get_Params()
        theta  = self.get_Theta()
        delta  = float(self.components.SlabDelta.text)
        layers = self.components.LayerList.items
        layers = self.sort_LayerList(layers,reverse=False)
        #
        mlayers = []
        for layer in layers:
            comp = compound.parse_species_formula(str(layer[1]))
            if comp == []:
                print "Error parsing formula"
                return
            thick   = float(layer[2])
            density = float(layer[3])
            rough   = float(layer[4])
            m = interface_model.Layer(comp=comp,density=density,
                                     roughness=rough,thickness=thick)
            mlayers.append(copy.copy(m))
        #
        model = interface_model.Model(substrate=None,
                                     layers=mlayers,
                                     top=None,
                                     theta=theta,
                                     params=params)
        model.slabify(delta=delta)
        model.set_param(fyel=params['fyel'])
        #
        if model:
            self.set_data(mname,model)
            self.components.Model.text = mname

    ######################################################

    ###########################################################
    #   Group, Model, Data Components etc
    ###########################################################
    
    ######################################################

    def init_Grp(self):
        " Initialize the menu    "
        grp = self.components.Grp.text
        tmp = self.shell.interp.symbol_table.list_symbols(tunnel=False)
        self.components.Grp.items = tmp['ins']
        self.components.Grp.text = grp
        return

    def on_Grp_select(self, event):
        "Re-init Node list given the new grp name"
        grp = self.components.Grp.stringSelection
        self.components.Grp.text = grp
        self.init_Model()
        return

    def on_Grp_keyDown(self,event):
        "select a variable name and check it out"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_Model()
        else:
            event.skip()
        return

    def init_Model(self):
        """
        Initialize the menu. Use the group thats
        selected in the group menu
        """
        grp = self.components.Grp.text  
        if len(grp) == 0: grp = None
        model = self.components.Model.text
        tmp = self.shell.interp.symbol_table.list_symbols(symbol=grp,tunnel=False)
        tmp = tmp['var'] + tmp['ins']
        tmp.sort()
        self.components.Model.items = tmp
        if model in tmp:
            self.components.Model.text = model
        else:
            self.components.Model.text = ''
        return

    def on_Model_select(self,event):
        "select a variable name and check it out"
        self.init_GUI()
        return

    def on_Model_keyDown(self,event):
        "select a variable name and check it out"
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_GUI()
        else:
            event.skip()
        return

    def get_ModelName(self,):
        if len(self.components.Model.stringSelection) > 0:
            self.components.Model.text = self.components.Model.stringSelection
        model = self.components.Model.text
        if len(model.strip()) == 0: return None
        name = "%s" % model.strip()
        return name

    def get_Model(self):
        name = self.get_ModelName()
        if name == None: return None
        return self.get_data(name)

    def set_Model(self,model):
        name = self.get_ModelName()
        return self.set_data(name,model)

    def CheckModelVar(self,):
        try:
            name = self.get_ModelName()
            m    = self.get_data(name)
            if type(m) == types.InstanceType:
                if hasattr(m,'layer'):
                    self.post_message("Valid model object: %s" % name)
                    self.components.SaveName.text = name
                    return True
                else:
                    self.post_message("Invalid model object: %s" % name)
                    return False
            else:
                self.post_message("Invalid model object")
                return False
        except:
            self.post_message("Invalid model object")
            return False

    ######################################################

    ###########################################################
    #   Layer parameters
    ###########################################################

    ######################################################

    def init_LayerParams(self,):
        self.components.LayerEdit.text   = ''
        self.components.CompEdit.text    = ''
        self.components.ThickEdit.text   = ''
        self.components.DensityEdit.text = ''
        self.components.RoughEdit.text   = ''
        
    def put_LayerParamFields(self,params):
        self.components.LayerEdit.text   = params[0]
        self.components.CompEdit.text    = params[1]
        self.components.ThickEdit.text   = params[2]
        self.components.DensityEdit.text = params[3]
        self.components.RoughEdit.text   = params[4]

    #def on_LayerEdit_keyDown(self,event):
    #    keyCode = event.keyCode
    #    if keyCode == wx.WXK_RETURN:
    #        self.update_LayerParams()
    #    else:
    #        event.skip()
    #    return

    def on_CompEdit_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_LayerParams()
        else:
            event.skip()
        return

    def on_ThickEdit_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_LayerParams()
        else:
            event.skip()
        return

    def on_DensityEdit_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_LayerParams()
        else:
            event.skip()
        return

    def on_RoughEdit_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_LayerParams()
        else:
            event.skip()
        return

    def update_LayerParams(self):
        param = ['','','','','']
        param[0] = self.components.LayerEdit.text.strip()
        if len(param[0]) == 0: return
        param[1] = self.components.CompEdit.text.strip()
        param[2] = self.components.ThickEdit.text.strip()
        param[3] = self.components.DensityEdit.text.strip()
        param[4] = self.components.RoughEdit.text.strip()
        #
        tmp = self.components.LayerList.items
        tmp = self.sort_LayerList(tmp,reverse=False)
        #
        found = False
        for j in range(len(tmp)):
            if tmp[j][0] == param[0]:
                tmp[j] = param
                found = True
                break
        if found == False:
            if param[0] == '0':
                tmp.insert(0,param)
            else:
                tmp.append(param)
            for j in range(len(tmp)):
                tmp[j][0] = str(j)
        tmp = self.sort_LayerList(tmp,reverse=True)
        self.components.LayerList.items = tmp

    ######################################################

    ###########################################################
    #   Layer List
    ###########################################################

    ######################################################
    def init_LayerList(self,):
        self.components.LayerList.items = []
    
    def on_LayerList_select(self,event):
        # print self.components.PkParams.items
        selected = self.components.LayerList.getStringSelection()
        self.put_LayerParamFields(selected[0])
        return
            
    def on_LayerInsert_mouseClick(self,event):
        layer = ['{(1.0)[N2]}','1000.','0.001','0.0']
        tmp  = self.components.LayerList.items
        nlayer = len(tmp)
        layer.insert(0,str(nlayer))
        tmp.append(layer)
        tmp = self.sort_LayerList(tmp)
        self.components.LayerList.items = tmp
        self.components.LayerList.SetSelection(0)

    def on_LayerDelete_mouseClick(self,event):
        selected = self.components.LayerList.getStringSelection()
        if selected == []: return
        sidx = int(selected[0][0])
        
        llist  = self.components.LayerList.items
        sorted = self.sort_LayerList(llist,reverse=False)
        
        sorted.pop(sidx)
        for j in range(len(sorted)):
            sorted[j][0] = str(j)

        sorted = self.sort_LayerList(sorted)
        self.components.LayerList.items = sorted

        if len(sorted) == 0:
            self.init_LayerParams()
        else:
            sel = (len(sorted)-1) - (sidx)
            if sel < 0: sel = 0
            self.components.LayerList.SetSelection(sel)
      
    def on_LayerMoveUp_mouseClick(self,event):
        selected = self.components.LayerList.getStringSelection()
        if selected == []: return
        selected = selected[0]
        sidx = int(selected[0])
        
        llist  = self.components.LayerList.items
        sorted = self.sort_LayerList(llist,reverse=False)
        if sidx == len(llist)-1: return
        
        sorted.pop(sidx)
        sorted.insert(sidx+1,selected)
        for j in range(len(sorted)):
            sorted[j][0] = str(j)

        sorted = self.sort_LayerList(sorted)
        self.components.LayerList.items = sorted
        
        sel = (len(sorted)-1) - (sidx+1)
        self.components.LayerList.SetSelection(sel)

    def on_LayerMoveDown_mouseClick(self,event):
        selected = self.components.LayerList.getStringSelection()
        if selected == []: return
        selected = selected[0]
        sidx = int(selected[0])
        
        llist  = self.components.LayerList.items
        sorted = self.sort_LayerList(llist,reverse=False)
        if sidx == 0: return
        
        sorted.pop(sidx)
        sorted.insert(sidx-1,selected)
        for j in range(len(sorted)):
            sorted[j][0] = str(j)

        sorted = self.sort_LayerList(sorted)
        self.components.LayerList.items = sorted

        sel = (len(sorted)-1) - (sidx-1)
        self.components.LayerList.SetSelection(sel)

    def sort_LayerList(self,layer,reverse=True):
        """
        given list re-order from top to bottom
        """
        xx  = {}
        layer = copy.copy(layer)
        for ll in layer:
            lidx = ll.pop(0)
            xx.update({lidx:ll})
        idx = range(len(xx))
        #
        if reverse==True: idx.reverse()
        zz = []
        for j in idx:
            yy = [str(j)] + xx[str(j)]
            zz.append(yy)
        return zz

    ######################################################

    ###########################################################
    #   Data
    ###########################################################

    ######################################################
    def init_Data(self):
        self.init_Theta()
        self.init_R_FY()

    ######################################################
    def init_Theta(self):
        """
        init theta
        """
        grp = self.components.Grp.text
        model = self.components.Model.text
        if len(grp) == 0: grp = None
        thet = self.components.Theta.text
        #
        tmp = self.shell.interp.symbol_table.list_symbols(symbol=grp,tunnel=False)
        tmp = tmp['var']
        if len(model)>0:
            tmp2 = self.shell.interp.symbol_table.list_symbols(symbol=model,tunnel=False)
            tmp  = tmp + tmp2['var']
        tmp.sort()
        tmp.insert(0,'Num.arange(0.01,1.01,0.01)')
        self.components.Theta.items = tmp
        #
        if thet in tmp:
            self.components.Theta.text = thet
        else:
            #self.components.Theta.text = 'Num.arange(0.01,1.01,0.01)'
            self.components.Theta.text = ''
        return

    ######################################################
    def init_R_FY(self):
        """
        init R and FY
        """
        grp   = self.components.Grp.text
        model = self.components.Model.text
        if len(grp) == 0: grp = None
        R   = self.components.RData.text
        FY  = self.components.FYData.text
        tmp = self.shell.interp.symbol_table.list_symbols(symbol=grp,tunnel=False)
        tmp = tmp['var']
        if len(model)>0:
            tmp2 = self.shell.interp.symbol_table.list_symbols(symbol=model,tunnel=False)
            tmp  = tmp + tmp2['var']
        tmp.sort()
        tmp.insert(0,'')
        self.components.RData.items = tmp
        self.components.FYData.items = tmp
        if R in tmp:
            self.components.RData.text = R
        else:
            self.components.RData.text = ''
        if FY in tmp:
            self.components.FYData.text = FY
        else:
            self.components.FYData.text = ''
        return

    ######################################################
    def get_Theta(self):
        thetstr = self.components.Theta.text
        thetstr = thetstr.strip()
        if len(thetstr) == 0: return None
        theta = self.eval_line(thetstr)
        return theta

    ######################################################
    def get_Rdat(self):
        Rstr = self.components.RData.text
        Rstr = Rstr.strip()
        if len(Rstr) == 0: return None
        R = self.eval_line(Rstr)
        return R

    ######################################################
    def get_FYdat(self):
        Ystr = self.components.FYData.text
        Ystr = Ystr.strip()
        if len(Ystr) == 0: return None
        Y = self.eval_line(Ystr)
        return Y

    ######################################################
    def on_DataUpdate_mouseClick(self,event):
        model = self.get_Model()
        if model == None: return

        theta = self.get_Theta()
        model.set_theta(theta)
        self.UpdateGuiFromModel()
        
        #ref  = self.get_Rdat()
        #fy  = self.get_FYdat()
        #...

    ######################################################

    ###########################################################
    #   Parameters
    ###########################################################

    ######################################################
    def init_Params(self):
        self.components.Energy.text = '10000.0'
        self.components.ConvWidth.text = '0.02'
        self.components.SampLen.text = '50.0'
        self.components.BeamVert.text = '0.01'
        self.components.BeamHorz.text = '10.0'
        self.components.AreaFlag.stringSelection = '0.0'
        self.components.RefScale.text = '1.0'
        self.components.FyEl.text = '0.0'
        self.components.FyEnergy.text = '6000.0'
        self.components.DetAngle.text = '90.0'
        self.components.ThetaNorm.text = '1.0'
        self.components.RoughFlag.stringSelection = '1.0'
        self.components.DelZ.text  = '10.0'
        self.components.PDepth.text = '3.0'
    
    ######################################################
    def get_Params(self):
        params = {}
        params['energy']   = float(self.components.Energy.text)
        params['wconv']    = float(self.components.ConvWidth.text)
        params['slen']     = float(self.components.SampLen.text)
        params['bvert']    = float(self.components.BeamVert.text)
        params['bhorz']    = float(self.components.BeamHorz.text)
        params['aflag']    = float(self.components.AreaFlag.stringSelection)
        params['rscale']   = float(self.components.RefScale.text)
        try:
            params['fyel']    = float(self.components.FyEl.text)
        except:
            params['fyel']    = str(self.components.FyEl.text)
        try:
            params['fyenergy'] = float(self.components.FyEnergy.text)
        except:
            params['fyenergy'] = str(self.components.FyEnergy.text)
        #
        params['adet']     = float(self.components.DetAngle.text)
        params['tnorm']    = float(self.components.ThetaNorm.text)
        params['rflag']    = float(self.components.RoughFlag.stringSelection)
        params['delz']     = float(self.components.DelZ.text)
        params['pdepth']   = float(self.components.PDepth.text)
        #
        return params

    ######################################################

    ###########################################################
    #   Calc/plot
    ###########################################################

    ######################################################
    def on_CalcR_mouseClick(self,event):
        self.UpdateModelFromGui()
        self.calc_R()
        
    def calc_R(self):
        model = self.get_Model()
        if model == None:
            print "Please save the model before calculating"
            return
        #
        if self.components.ShowTime.checked:
            t = time.time()
            model.calc_R()
            print " Elapsed time to calc R = %f seconds\n" % (time.time() - t)
        else:
            model.calc_R()
        #
        self.AutoPlot(FYflag=False)

    ######################################################
    def on_CalcFY_mouseClick(self,event):
        self.UpdateModelFromGui()
        self.calc_FY()

    def calc_FY(self):
        model = self.get_Model()
        if model == None:
            print "Please save the model before calculating"
            return
        #
        if self.components.ShowTime.checked:
            t = time.time()
            model.calc_FY()
            print " Elapsed time to calc FY = %f seconds\n" % (time.time() - t)
        else:
            model.calc_FY()
        #
        self.AutoPlot()

    #############################################################
    def AutoPlot(self,FYflag=True):
        model = self.get_Model()
        if model == None: return
        #
        hold = self.components.HoldRFY.checked
        bar  = self.components.BarPlot.checked
        if self.components.PlotRFY.checked:
            if self.components.PlotData.checked:
                Rdat  = self.get_Rdat()
                FYdat = self.get_FYdat()
            else:
                Rdat = None
                FYdat = None
            model.plot(ty='calc',hold=hold,FYflag=FYflag,
                       Rdat=Rdat,FYdat=FYdat)
        if self.components.PlotDensity.checked:
            model.plot(ty='density',bar=bar,hold=hold)
        if self.components.PlotComps.checked:
            model.plot(ty='comp',bar=bar,hold=hold)
        if self.components.PlotElements.checked:
            model.plot(ty='el',bar=bar,hold=hold)
        if self.components.PlotFracs.checked:
            model.plot(ty='frac',bar=bar,hold=hold)

##################################################
if __name__ == '__main__':
    app = model.Application(wxXrrBuilder)
    app.MainLoop()
