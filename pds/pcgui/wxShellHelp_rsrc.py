{'application':{'type':'Application',
          'name':'Template',
    'backgrounds': [
    {'type':'Background',
          'name':'bgTemplate',
          'title':u'PDS Shell Help',
          'size':(638, 724),
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
         ]
     },
         'components': [

{'type':'TextArea', 
    'name':'HelpText', 
    'position':(10, 10), 
    'size':(612, 643), 
    'editable':False, 
    },

] # end components
} # end background
] # end backgrounds
} }
