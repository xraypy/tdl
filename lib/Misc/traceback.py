import sys
import traceback
from Util import PrintExceptErr, PrintShortExcept, SymbolError, ConstantError

def testx():
    x = 1
    if x < 10:
        raise SymbolError, 'oops'
    return 'hello'


try:
    testx()
except SymbolError, detail:
    print 'bad:: ', detail
    traceback.print_exc()
    


