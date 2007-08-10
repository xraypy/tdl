###################################################################
def testCommand2Expr():
    import Symbol
    s  = Symbol.SymbolTable()
    s.addVariable('avar',3)
    s.addVariable('bvar',3)
    s.addVariable('x',Num.arange(10))
    print Command2Expr('func avar, bvar cvar, dvar=3',symtable=s)
    print Command2Expr('func *.*  cvar, dvar=3',symtable=s)
    print Command2Expr('plot x x*2',symtable=s)
    
