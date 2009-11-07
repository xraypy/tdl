main = """This is tdl main help"""

import tdlHelpTopics
Help_topics = tdlHelpTopics.generate()

class Helper(object):
    """Helper looks up an displays help topics
    and/or pydoc help on tdl/python objects"""
    
    def __init__(self,*args,**kws):
        self.tdl = None
        self.buff = []

    def help(self,*args,**kws):
        "return help text for a series of arguments"
        for arg in args:
            if arg is None: continue
            # print arg, arg in Help_topics
            if arg in Help_topics:
                # print Help_topics[arg]
                self.addtext(Help_topics[arg])
            else:
                sym = self.tdl.symtable.getSymbol(arg,create=False)
                self.addtext("   %s = %s" % (arg, repr(sym)))
                
                
    def addtext(self,text):
        self.buff.append(text)

    def getbuffer(self,delim='\n'):
        out = delim.join(self.buff)
        self.buff = []
        return out
        

#     "show help on topic or object"
#     outbuff = []
#     has_tdl = tdl is not None
# 
#     for a in args:
#         
#             outbuff.append(pydoc.help(a))
#     else:
#         
#     try:
#         f = open(name)
#         l = f.readlines()
# 
#     except IOError:
#         print "cannot open file: %s." % name
#         return
#     finally:
#         f.close()
#     show_more(l,filename=name,pagelength=pagelength)
