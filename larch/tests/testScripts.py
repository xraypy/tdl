try:
    import lib as larch
except:
    import larch
import sys
import os

interp = larch.interp()
input  = larch.input()

def doFile(inputfile,output=None):
    try:
        input.readfile(inputfile)
    except IOError:
        print "could not read %s" % inputfile

    def out(s):  print s
    save_stdout = sys.stdout
    save_stdout.write( 'testing %s ..' % inputfile)

    if output is not None:
        fout = open(output,'w')
        sys.stdout = fout
        def out(s):
            fout.write("%s\n"% s)
        
    while input:
        block,fname,lineno = input.get()
        ret = interp.eval(block,fname=fname,lineno=lineno)
        if callable(ret) and not isinstance(ret,type):
            try:
                if 1 == len(block.split()):
                    ret = ret()
            except:
                pass
        if interp.error:
            out('== Error ==')
            err  = interp.error.pop(0)
            out( "%s: %s" % err.get_error())
            for err in interp.error:
                err_type,err_msg =  err.get_error()
                if not (err_type.startswith('Extra Error') or
                        err_type.startswith('Eval Error') ):
                    out( err_msg)
            out('===========')
        if ret is not None:
            out(ret)

    fout.flush()
    fout.close()
    sys.stdout = save_stdout 
    
def compareFiles(fname1,fname2):
    
    text1 = open(fname1,'r').readlines()
    text2 = open(fname2,'r').readlines()

    diff = []
    if len(text1) != len(text2):
        diff.append("files have different lengths: %i, %i " % (len(text1),
                                                               len(text2)))
    i = 0
    for t1,t2 in zip(text1,text2):
        i = i+1
        if t1 != t2:
            diff.append("line %i:\n expected '%s'\n got    's'" % (t1,t2))
    if len(diff)>0:
        diff.insert(0,"%s and %s differ!"%(fname1,fname2))
    if len(diff)>0:
        return "\n".join(diff)
    return ' OK.' 
            

def testScript(script,output):
    tmpfile = '%s.test'% script
    doFile(script,output=tmpfile)
    print compareFiles(output, tmpfile)
    os.unlink(tmpfile)
    
if __name__ == '__main__':
    for t in ('t1','t2'):
        script = '%s.lar' % t
        output = '%s.out' % t
        testScript(script,output)
