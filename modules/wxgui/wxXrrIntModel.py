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
import numpy as num

from   wxUtil import wxUtil
import interface_model

#########################################################################

class wxXrrIntModel(model.Background, wxUtil):

    ###########################################################
    # Init and util methods
    ###########################################################
    
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.components.DistList._autoresize = 0
        self.startup  = True
        self.dir      = '.'

        # set up shell
        self.shell = None
        self.init_shell()

        # Make sure interface_model is loaded
        self.exec_line("import interface_model")

        # init all the menus
        self.init_ShellItems()
        self.init_GUI()
        self.current = True
        
   ###########################################################

    def init_ShellItems(self,):
        self.init_Grp()
        self.init_Model()
        self.init_Data()

    def on_UpdateShell_mouseClick(self,event):
        self.init_ShellItems()

    def on_Clear_mouseClick(self,event):
        self.components.Grp.text = ""
        self.components.Model.text = ""
        self.init_ShellItems()
        self.init_GUI()
    
    def init_GUI(self):
        if self.startup:
            self.startup = False
        self.init_Components()
        self.init_DistList()
        self.init_DistParams()
        self.init_Params()

        check = self.CheckModelVar()
        if check == False:
            return
        else:
            self.UpdateGuiFromModel()

    ###########################################################
    #   Menu
    ###########################################################

    def on_menuFileExit_select(self,event): 
        self.close()

    def on_menuHelpParams_select(self,event): 
        import wxXrrIntModelHelp
        wxXrrIntModelHelp = mod_import(wxXrrIntModelHelp)
        dir       = os.path.dirname(wxXrrIntModelHelp.__file__)
        filename  = os.path.join(dir,'wxXrrIntModelHelp.rsrc.py')
        wxXrrIntModelHelp = wxXrrIntModelHelp.wxXrrIntModelHelp
        self.wxXrrIntModelHelp = model.childWindow(self,wxXrrIntModelHelp,
                                                   filename=filename)
        self.wxXrrIntModelHelp.position = (200, 5)
        self.wxXrrIntModelHelp.visible = True

    ######################################################

    ###########################################################
    #   Update GUI / Update Model
    ###########################################################

    ######################################################
    def UpdateGuiFromModel(self,):
        """
        update gui from a model
        """
        model = self.get_Model()
        if model == None: return

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
            #
            if params.has_key('fyel'):
                self.components.FyEl.text   = str(params['fyel'])
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

        # components
        self.init_Components()

        # Dist params        
        self.init_DistParamsFromModel()

    ######################################################
    def UpdateModelDistParamsFromGui(self):
        """
        par[0] = index
        par[1] = type
        par[2] = zst
        par[3] = zen
        par[4] = CX (param1)
          :
        """
        dlist = copy.copy(self.components.DistList.items)
        dlist = self.sort_DistList(dlist,reverse=False)
        inter = []
        try:
            for par in dlist:
                dist = {}
                type = par[1]
                dist['type'] = type
                #
                zst  = par[2]
                zen  = par[3]
                if len(zst)>0:
                    dist['zst'] = eval(zst)
                if len(zen)>0:
                    dist['zen'] = eval(zen)
                #
                CX   = float(par[4])
                dist['CX']  = CX
                #
                if dist['type'] in ('erf','erfc','exp','expc','gauss'):
                    cen  = par[5]
                    if (eval(cen) == None):
                        dist['cen'] = 0.0
                    elif (len(cen) == 0):
                        dist['cen'] = 0.0
                    else:
                        dist['cen'] = float(cen)
                    #
                    sig  = par[6]
                    if (eval(sig) == None):
                        dist['sig'] = 0.0
                    elif (len(sig) == 0):
                        dist['sig'] = 0.0
                    else:
                        dist['sig'] = float(sig)
                #
                elif dist['type'] == 'linear':
                    CXen  = par[5]
                    if (eval(CXen) == None):
                        dist['CXen'] = 0.0
                    elif (len(CXen) == 0):
                        dist['CXen'] = 0.0
                    else:
                        dist['CXen'] = float(CXen)
                #
                inter.append(copy.copy(dist))
        except:
            print "Error setting distribution parameter"
            return None
        
        #########

        # norm (just for selected comp)...
        norm  = int(self.components.Norm.stringSelection)
        scale = self.components.ScaleToNorm.checked

        # density (for all distr)...
        rhoflag = int(self.components.DensityFlag.stringSelection)
        rhoscale = self.components.DensityScale.checked

        #########
        model = self.get_Model()
        if model == None: return
        comp = str(self.components.Components.stringSelection)
        if len(comp.strip()) == 0:
            print "Error getting component dist params"
            return
        idx = model.slab._comp_idx(comp)
        if idx == -1:
            print "Error getting component dist params"
            return
        model.slab.distpar[idx].inter = inter
        model.slab.distpar[idx].norm  = norm
        model.slab.distpar[idx].scale_to_norm = scale
        #
        model.slab.rhoflag = rhoflag
        model.slab.rhoscale = rhoscale

        ########
        self.init_DistParamsFromModel()
        if self.components.AutoCalcDist.checked:
            self.calc_dist()

    ######################################################

    ###########################################################
    #   Group, Model Name, etc
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
                if hasattr(m,'ref'):
                    self.post_message("Valid model object: %s" % name)
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
    #   Components / density constraint
    ###########################################################

    ######################################################

    def init_Components(self):
        tmp = self.components.Components.stringSelection
        model = self.get_Model()
        if model == None:
            self.components.Components.items = []
            return
        # components
        comp = model.slab.list_comps()
        self.components.Components.items = comp
        if tmp in comp:
            self.components.Components.stringSelection = tmp
        else:
            self.components.Components.stringSelection = comp[0]

    def on_Norm_select(self,event):
        self.UpdateModelDistParamsFromGui()
    
    def on_ScaleToNorm_mouseClick(self,event):
        self.UpdateModelDistParamsFromGui()

    def on_Components_select(self,event):
        self.init_DistParamsFromModel()
        self.init_ParamSliders()

    def on_DensityFlag_select(self,event):
        self.UpdateModelDistParamsFromGui()
    
    def on_DensityScale_mouseClick(self,event):
        self.UpdateModelDistParamsFromGui()


    ######################################################

    ###########################################################
    #   Dist parameters 
    ###########################################################

    ######################################################

    def init_DistParams(self,):
        self.components.DistIndex.text = ''
        self.components.DistModel.stringSelection = 'box'
        self.components.Zstart.text = '' 
        self.components.Zend.text   = ''
        self.components.Param1.text = ''
        self.components.Param2.text = ''
        self.components.Param3.text = ''
        self.components.Param4.text = ''
        self.components.Param5.text = ''
        #
        self.components.param1Label.text = 'Param 1'
        self.components.param2Label.text = 'Param 2'
        self.components.param3Label.text = 'Param 3'
        self.components.param4Label.text = 'Param 4'
        self.components.param5Label.text = 'Param 5'
        #
        self.init_ParamSliders()
        
    ##################################################
    def init_ParamSliders(self):
        self.components.P1slider.value = 0.
        self.components.P2slider.value = 0.
        self.components.P3slider.value = 0.
        self.components.P4slider.value = 0.
        self.components.P5slider.value = 0.
        #
        self.p1 = self.components.Param1.text
        self.p2 = self.components.Param2.text
        self.p3 = self.components.Param3.text
        self.p4 = self.components.Param4.text
        self.p5 = self.components.Param5.text
        #
        self.p1cnt = 0
        self.p2cnt = 0
        self.p3cnt = 0
        self.p4cnt = 0
        self.p5cnt = 0
        #
        
    ##################################################
    def put_DistParamFields(self,params):
        self.components.DistIndex.text = params[0]
        self.components.DistModel.stringSelection = params[1]
        self.components.Zstart.text = params[2] 
        self.components.Zend.text   = params[3]
        self.components.Param1.text = params[4]
        self.components.Param2.text = params[5]
        self.components.Param3.text = params[6]
        self.components.Param4.text = params[7]
        self.components.Param5.text = params[8]
        #
        type = self.components.DistModel.stringSelection
        if type == 'box':
            self.components.param1Label.text = 'CX'
            self.components.param2Label.text = 'None'
            self.components.param3Label.text = 'None'
            self.components.param4Label.text = 'None'
            self.components.param5Label.text = 'None'
        elif type in ('erf','erfc','exp','expc','gauss'):
            self.components.param1Label.text = 'CX'
            self.components.param2Label.text = 'Center'
            self.components.param3Label.text = 'Width'
            self.components.param4Label.text = 'None'
            self.components.param5Label.text = 'None'
        elif type == 'linear':
            self.components.param1Label.text = 'CX'
            self.components.param2Label.text = 'CXen'
            self.components.param3Label.text = 'None'
            self.components.param4Label.text = 'None'
            self.components.param5Label.text = 'None'
        else:
            self.components.param1Label.text = 'Param 1'
            self.components.param2Label.text = 'Param 2'
            self.components.param3Label.text = 'Param 3'
            self.components.param4Label.text = 'Param 4'
            self.components.param5Label.text = 'Param 5'

    ##################################################
    def on_DistModel_select(self,event):
        self.update_DistParams()

    def on_Zstart_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_DistParams()
        else:
            event.skip()
        return

    def on_Zend_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_DistParams()
        else:
            event.skip()
        return

    ##################################################
    def on_Param1_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_DistParams()
            self.init_ParamSliders()
        else:
            event.skip()
        return

    def on_Param2_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_DistParams()
            self.init_ParamSliders()
        else:
            event.skip()
        return

    def on_Param3_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_DistParams()
            self.init_ParamSliders()
        else:
            event.skip()
        return

    def on_Param4_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_DistParams()
            self.init_ParamSliders()
        else:
            event.skip()
        return

    def on_Param5_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.update_DistParams()
            self.init_ParamSliders()
        else:
            event.skip()
        return

    def on_SliderScale_keyDown(self,event):
        keyCode = event.keyCode
        if keyCode == wx.WXK_RETURN:
            self.init_ParamSliders()
        else:
            event.skip()
        return

    ##################################################
    def on_P1slider_select(self,event):
        if self.p1cnt < 2:
            self.p1cnt = self.p1cnt + 1
            return
        else:
            self.p1cnt = 1
        try:
            #delta = self.components.P1slider.value
            delta = float(event.target.value)
            scale = float(self.components.SliderScale.text)
            old   = float(self.p1)
        except:
            self.init_ParamSliders()
            return
        new   = old + old*(scale/100.) * (delta/100.) 
        new   = str(new)
        self.components.Param1.text = new
        self.update_DistParams()

    def on_P2slider_select(self,event):
        if self.p2cnt < 2:
            self.p2cnt = self.p2cnt + 1
            return
        else:
            self.p2cnt = 1
        try:
            #delta = self.components.P1slider.value
            delta = float(event.target.value)
            scale = float(self.components.SliderScale.text)
            old   = float(self.p2)
        except:
            self.init_ParamSliders()
            return
        new   = old + old*(scale/100.) * (delta/100.) 
        new   = str(new)
        self.components.Param2.text = new
        self.update_DistParams()

    def on_P3slider_select(self,event):
        if self.p3cnt < 2:
            self.p3cnt = self.p3cnt + 1
            return
        else:
            self.p3cnt = 1
        try:
            #delta = self.components.P1slider.value
            delta = float(event.target.value)
            scale = float(self.components.SliderScale.text)
            old   = float(self.p3)
        except:
            self.init_ParamSliders()
            return
        new   = old + old*(scale/100.) * (delta/100.) 
        new   = str(new)
        self.components.Param3.text = new
        self.update_DistParams()

    def on_P4slider_select(self,event):
        if self.p4cnt < 2:
            self.p4cnt = self.p4cnt + 1
            return
        else:
            self.p4cnt = 1
        try:
            #delta = self.components.P1slider.value
            delta = float(event.target.value)
            scale = float(self.components.SliderScale.text)
            old   = float(self.p4)
        except:
            self.init_ParamSliders()
            return
        new   = old + old*(scale/100.) * (delta/100.) 
        new   = str(new)
        self.components.Param4.text = new
        self.update_DistParams()

    def on_P5slider_select(self,event):
        if self.p5cnt < 2:
            self.p5cnt = self.p5cnt + 1
            return
        else:
            self.p5cnt = 1
        try:
            #delta = self.components.P1slider.value
            delta = float(event.target.value)
            scale = float(self.components.SliderScale.text)
            old   = float(self.p5)
        except:
            self.init_ParamSliders()
            return
        new   = old + old*(scale/100.) * (delta/100.) 
        new   = str(new)
        self.components.Param5.text = new
        self.update_DistParams()

    ##################################################
    def update_DistParams(self):
        """
        update the list of dist params
        """
        # make sure a component is selected
        comp = self.components.Components.stringSelection
        if len(comp.strip()) == 0: return
        #
        param = ['','','','','','','','','']
        idx      = self.components.DistIndex.text.strip()
        param[0] = idx
        if len(param[0]) == 0: return
        param[1] = self.components.DistModel.stringSelection.strip()
        param[2] = self.components.Zstart.text.strip()
        param[3] = self.components.Zend.text.strip()
        param[4] = self.components.Param1.text.strip()
        param[5] = self.components.Param2.text.strip()
        param[6] = self.components.Param3.text.strip()
        param[7] = self.components.Param4.text.strip()
        param[8] = self.components.Param5.text.strip()
        #
        tmp = self.components.DistList.items
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
        #
        tmp = self.sort_DistList(tmp,reverse=True)
        self.components.DistList.items = tmp
        self.UpdateModelDistParamsFromGui()

    ######################################################
    def init_DistParamsFromModel(self):
        #
        model = self.get_Model()
        if model == None: return
        comp = str(self.components.Components.stringSelection)
        if len(comp.strip()) == 0:
            print "Error getting component dist params"
            return
        idx = model.slab._comp_idx(comp)
        if idx == -1:
            print "Error getting component dist params"
            return
        dpar = model.slab.distpar[idx]

        # density
        rhoflag = str(int(model.slab.rhoflag))
        self.components.DensityFlag.stringSelection = rhoflag
        self.components.DensityScale.checked = model.slab.rhoscale

        # norm
        norm = str(int(dpar.norm))
        self.components.Norm.stringSelection = norm
        
        # scale to norm
        scale = dpar.scale_to_norm
        self.components.ScaleToNorm.checked = scale

        #
        top = "type=%s, zst=%6.3f, zen =%6.3f, CX=%g" %(dpar.top['type'],
                                                        dpar.top['zst'],
                                                        dpar.top['zen'],
                                                        dpar.top['CX'])
        self.components.Top.text = top
        subs = "type=%s, zst=%6.3f, zen =%6.3f, CX=%g" %(dpar.subs['type'],
                                                        dpar.subs['zst'],
                                                        dpar.subs['zen'],
                                                        dpar.subs['CX'])
        self.components.Subs.text = subs
        
        #
        dlist = []
        for j in range(len(dpar.inter)):
            param = []
            param.append(str(j))
            param.append(str(dpar.inter[j].get('type')))
            param.append(str(dpar.inter[j].get('zst')))
            param.append(str(dpar.inter[j].get('zen')))
            param.append(str(dpar.inter[j].get('CX')))
            # below may depend on type...
            if param[1] == 'box':
                param.append(str(None))
                param.append(str(None))
                param.append(str(None))
                param.append(str(None))
            elif param[1] in ('erf','erfc','exp','expc','gauss'):
                param.append(str(dpar.inter[j].get('cen')))
                param.append(str(dpar.inter[j].get('sig')))
                param.append(str(None))
                param.append(str(None))
            elif param[1] == 'linear':
                param.append(str(dpar.inter[j].get('CXen')))
                param.append(str(None))
                param.append(str(None))
                param.append(str(None))
            #
            dlist.append(copy.copy(param))
        #######
        # post dlist and set selection
        #######
        sidx = self.components.DistIndex.text.strip()
        npar = len(dlist)
        dlist = self.sort_DistList(dlist,reverse=True)
        self.components.DistList.items = dlist
        try:
            sidx = npar - (int(sidx) + 1)
        except:
            sidx = 0
        if (sidx >= 0) and (sidx < len(dlist)-1):
            self.components.DistList.SetSelection(sidx)
        else:
            self.components.DistList.SetSelection(0)
    
    ##################################################

    ###########################################################
    #   Dist List
    ######################################################

    ##################################################
    def init_DistList(self,):
        self.components.Top.text = ''
        self.components.Subs.text = ''
        self.components.DistList.items = []
    
    def on_DistList_select(self,event):
        oldidx = self.components.DistIndex.text
        selected = self.components.DistList.getStringSelection()
        self.put_DistParamFields(selected[0])
        newidx = self.components.DistIndex.text
        if oldidx != newidx:
            self.init_ParamSliders()
        
    def on_DistInsert_mouseClick(self,event):
        # make sure a component is selected
        comp = self.components.Components.stringSelection
        if len(comp.strip()) == 0: return

        dist = ['box','None','None','0.0','0.0','0.0','0.0','0.0']
        tmp  = self.components.DistList.items
        ndist = len(tmp)
        dist.insert(0,str(ndist))
        tmp.append(dist)
        tmp = self.sort_DistList(tmp,reverse=True)
        self.components.DistList.items = tmp
        self.components.DistList.SetSelection(0)
        #
        self.UpdateModelDistParamsFromGui()

    def on_DistDelete_mouseClick(self,event):
        selected = self.components.DistList.getStringSelection()
        if selected == []: return
        sidx = int(selected[0][0])
        
        dlist  = self.components.DistList.items
        sorted = self.sort_DistList(dlist,reverse=False)
        
        sorted.pop(sidx)
        for j in range(len(sorted)):
            sorted[j][0] = str(j)

        self.components.DistList.items = sorted

        if len(sorted) == 0:
            self.init_DistParams()
        else:
            sel = sidx - 1
            if sel < 0: sel = 0
            self.components.DistList.SetSelection(sel)
        #
        self.UpdateModelDistParamsFromGui()
    
    def sort_DistList(self,dist,reverse=True):
        """
        given list re-order by first col
        """
        xx  = {}
        dist = copy.copy(dist)
        for dd in dist:
            didx = dd.pop(0)
            xx.update({didx:dd})
        idx = range(len(xx))
        if reverse==True: idx.reverse()
        #
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
        tmp.insert(0,'num.arange(0.01,1.01,0.01)')
        self.components.Theta.items = tmp
        #
        if thet in tmp:
            self.components.Theta.text = thet
        else:
            #self.components.Theta.text = 'num.arange(0.01,1.01,0.01)'
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
        params['adet']     = float(self.components.DetAngle.text)
        params['tnorm']    = float(self.components.ThetaNorm.text)
        params['rflag']    = float(self.components.RoughFlag.stringSelection)
        params['delz']     = float(self.components.DelZ.text)
        params['pdepth']   = float(self.components.PDepth.text)
        #
        return params

    ######################################################
    def on_ParamsUpdate_mouseClick(self,event):
        model = self.get_Model()
        if model == None: return
        #
        params = self.get_Params()
        model.set_param(**params)
        self.UpdateGuiFromModel()

    ######################################################

    ###########################################################
    #   Calc/plot
    ###########################################################

    ######################################################
    def on_CalcDist_mouseClick(self,event):
        self.calc_dist()
        
    def calc_dist(self):
        model = self.get_Model()
        if model == None: return
        #
        if self.components.ShowTime.checked:
            t = time.time()
            model.calc_dist()
            print " Elapsed time to calc dist = %f seconds\n" % (time.time() - t)
        else:
            model.calc_dist()
        #
        self.init_DistParamsFromModel()
        #
        if self.components.AutoCalcRFY.checked:
            if self.components.FyEl.text == '0.0':
                self.calc_R()
            else:
                self.calc_FY()
        else:
            self.AutoPlot()

    ######################################################
    def on_CalcR_mouseClick(self,event):
        self.calc_R()
        
    def calc_R(self):
        model = self.get_Model()
        if model == None: return
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
        self.calc_FY()

    def calc_FY(self):
        model = self.get_Model()
        if model == None: return
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
    
#############################################################################
#############################################################################
if __name__ == '__main__':
    app = model.Application(wxXrrModel)
    app.MainLoop()
