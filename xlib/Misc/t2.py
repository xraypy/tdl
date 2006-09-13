class Const:
    def __init__(self,val=None,desc=None):
        self.__x = val
    def _getx(self): return self.__x
    def _delx(self): self.__x = None
    
    def _setx(self,val):
        print ' SET X ', val
        if self.__x == None:
            self.__x = val
        else:
            print 'cannot set x!!'

    x = property(_setx,_getx,_delx,'')
