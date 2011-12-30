import numpy as num

# make some dat
def data1(en):
    pk1 = gauss(en,4.0,500.,0.5)
    pk2 = gauss(en,7.0,800.,0.6)
    bgr_pk = gauss(en,15.,100.,10.)
    data = 50.0 + .01*en + pk1 + pk2 + bgr_pk
    data = data + 20.0*num.random.rand(len(data))
    return(data)

def data2():
    data = 100
    return(data)

####
def gauss(energy,cen,ampl,fwhm):
    sigma = fwhm/2.35482
    return ampl*num.exp(-1*((cen-energy)**2.0) /(2.*sigma**2))
