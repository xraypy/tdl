data = {'application':{'type':'Application',
          'name':'Template',
    'backgrounds': [
    {'type':'Background',
          'name':'bgTemplate',
          'title':'Standard Template with File->Exit menu',
          'size':(526, 624),
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
    'size':(500, 536), 
    'editable':False, 
    },

] # end components
} # end background
] # end backgrounds
} }
