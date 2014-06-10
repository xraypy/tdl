data = {'application':{'type':'Application',
          'name':'Template',
    'backgrounds': [
    {'type':'Background',
          'name':'wxXRF_bgTemplate',
          'title':u'Scan Selector',
          'size':(586, 133),
          'style':['resizeable'],

        'menubar': {'type':'MenuBar',
         'menus': [
             {'type':'Menu',
             'name':'menuFile',
             'label':'&File',
             'items': [
                  {'type':'MenuItem',
                   'name':'menuFileExit',
                   'label':'E&xit',
                  },
              ]
             },
             {'type':'Menu',
             'name':'menuHelp',
             'label':u'Help',
             'items': [
                  {'type':'MenuItem',
                   'name':'menuHelpDocumentation',
                   'label':u'Documentation',
                  },
              ]
             },
         ]
     },
         'components': [

{'type':'Button', 
    'name':'Update', 
    'position':(405, 43), 
    'size':(72, 22), 
    'label':u'Update', 
    },

{'type':'StaticText', 
    'name':'StaticText1', 
    'position':(16, 16), 
    'text':u'Spec Path', 
    },

{'type':'TextField', 
    'name':'SpecPath', 
    'position':(76, 14), 
    'size':(347, -1), 
    },

{'type':'Button', 
    'name':'SpecPathSel', 
    'position':(429, 14), 
    'size':(48, 22), 
    'label':u'Select', 
    },

{'type':'StaticText', 
    'name':'Lbl111', 
    'position':(16, 48), 
    'text':u'Spec File(s)', 
    },

{'type':'Choice', 
    'name':'SpecFile', 
    'position':(76, 44), 
    'size':(273, -1), 
    'items':[], 
    },

{'type':'Button', 
    'name':'SpecFileSel', 
    'position':(354, 43), 
    'size':(48, 22), 
    'label':u'Select', 
    },

{'type':'Button', 
    'name':'Filter', 
    'position':(485, 23), 
    'size':(72, 35), 
    'label':u'Filter', 
    },

] # end components
} # end background
] # end backgrounds
} }
