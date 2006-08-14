{'application':{'type':'Application',
          'name':'Template',
    'backgrounds': [
    {'type':'Background',
          'name':'tdl_wxGUI_bgTemplate',
          'title':u'TDL',
          'size':(855, 916),
          'statusBar':1,
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
                   'command':'exit',
                  },
              ]
             },
             {'type':'Menu',
             'name':'menuShell',
             'label':'Shell',
             'items': [
                  {'type':'MenuItem',
                   'name':'menuShellStart',
                   'label':'Shell',
                  },
              ]
             },
             {'type':'Menu',
             'name':'menuWindow',
             'label':'Window',
             'items': [
                  {'type':'MenuItem',
                   'name':'menuWindowXRF',
                   'label':'XRF',
                  },
              ]
             },
         ]
     },
         'components': [

{'type':'StaticLine', 
    'name':'Sep1', 
    'position':(8, 770), 
    'size':(833, -1), 
    'font':{'faceName': 'Microsoft Sans Serif', 'family': 'sansSerif', 'size': 8}, 
    'layout':'horizontal', 
    },

{'type':'StaticText', 
    'name':'Prompt', 
    'position':(7, 787), 
    'font':{'faceName': 'Microsoft Sans Serif', 'family': 'sansSerif', 'size': 10}, 
    'text':u'>>>', 
    },

{'type':'TextField', 
    'name':'ShellCmd', 
    'position':(41, 783), 
    'size':(796, 27), 
    'font':{'faceName': 'Microsoft Sans Serif', 'family': 'sansSerif', 'size': 10}, 
    },

{'type':'TextArea', 
    'name':'ShellText', 
    'position':(9, 5), 
    'size':(829, 755), 
    'editable':False, 
    'font':{'faceName': 'Lucida Console', 'family': 'sansSerif', 'size': 9}, 
    },

] # end components
} # end background
] # end backgrounds
} }
