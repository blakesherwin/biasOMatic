import numpy as np
#from solenspipe import utility as simgen
#from solenspipe.utility import w_n
import pylab
#from falafel import utils
#import pytempura
import numpy as np
import pylab
import numpy
from scipy import stats
import sys
recompute = False#True
if sys.argv[1]:
    sample =sys.argv[1]
else:
    sample = 'red'
fgArray = ['0093_qe','0145_qe','0093_prh']#['0093_prh', '0093_psh', '0093_qe', '0145_prh', '0145_psh', '0145_qe', 'freqcoadd_prh', 'freqcoadd_psh', 'freqcoadd_qe', 'hilc-cibd_prh', 'hilc-cibd_psh', 'hilc-cibd_qe', 'hilc-tszd_prh', 'hilc-tszd_psh', 'hilc-tszd_qe', 'hilc_prh', 'hilc_psh', 'hilc_qe']
print(sample,fgArray)

#'freqcoadd_prh']#['0093_psh', '0145_prh', '0145_psh', 'freqcoadd_prh', 'freqcoadd_psh', 'freqcoadd_qe', 'hilc-cibd_prh', 'hilc-cibd_psh', 'hilc-cibd_qe', 'hilc-tszd_prh', 'hilc-tszd_psh', 'hilc-tszd_qe', 'hilc_prh', 'hilc_psh', 'hilc_qe']#['0145_qe','0093_qe','0093_prh']#, 'freqcoadd_prh', 'freqcoadd_psh', 'freqcoadd_qe', 'hilc-cibd_prh', 'hilc-cibd_psh', 'hilc-cibd_qe', 'hilc-tszd_prh', 'hilc-tszd_psh', 'hilc-tszd_qe', 'hilc_prh', 'hilc_psh', 'hilc_qe']###['0093_prh', '0093_psh', '0093_qe', '0145_prh', '0145_psh', '0145_qe']#, 'freqcoadd_prh', 'freqcoadd_psh', 'freqcoadd_qe', 'hilc-cibd_prh', 'hilc-cibd_psh', 'hilc-cibd_qe', 'hilc-tszd_prh', 'hilc-tszd_psh', 'hilc-tszd_qe', 'hilc_prh', 'h\
#ilc_psh', 'hilc_qe']        
outputDir="/global/homes/b/bsherwin/fgBiasOutputs/"

def bin_spectrum(ells, cells, bin_edges, ell_weighted=False):
    if ell_weighted:
        sums = stats.binned_statistic(ells, ells, statistic='sum', bins=bin_edges)
        cl = stats.binned_statistic(ells, ells * cells, statistic='sum', bins=bin_edges)
        cl = cl[0] / sums[0]
    else:
        cl = stats.binned_statistic(ells, cells, statistic='mean', bins=bin_edges)[0]

    return cl

ell_bin_edges = np.concatenate(([19.5], np.arange(51.5, 1000 + 50 / 2, 50)))
zmidarray =[]
clbiasarray =[]
import pickle

if recompute:
    from pixell import enmap, wcsutils,curvedsky as cs
    import healpy as hp
    mapCMBK = hp.read_map("/global/homes/b/bsherwin/kap.fits")
    clsRef = hp.anafast(mapCMBK,mapCMBK,lmax=1024)###cs.alm2cl(mapFG, mapK2)                                                                                                                                                                           
    ellsRef = numpy.arange(len(clsRef))
    clsRefBinned = bin_spectrum(ellsRef,clsRef,ell_bin_edges)/2.
    pickle.dump(clsRef,open(outputDir+'clsRef.pkl','wb'))
    pickle.dump(ellsRef,open(outputDir+'ellsRef.pkl','wb'))
    pickle.dump(clsRefBinned,open(outputDir+'clsRefBinned.pkl','wb'))
    print('done recomputing')
else:
    ellsRef = pickle.load(open(outputDir+'ellsRef.pkl','rb'))
    clsRef = pickle.load(open(outputDir+'clsRef.pkl','rb'))
    clsRefBinned = pickle.load(open(outputDir+'clsRefBinned.pkl','rb'))

b = numpy.zeros((23,len(clsRefBinned)))
br = numpy.zeros((23,len(clsRefBinned)))

for FGiter in range(len(fgArray)):
 count = 0
 FG = fgArray[FGiter]

 for i in numpy.arange(22.):
  if recompute:
    zone = i*0.2+0.2
    ztwo = zone + 0.1
    print(zone,ztwo)
    print("%.1f%.1f"%(zone,ztwo))
    mapK = hp.fitsfunc.read_map("/global/homes/b/bsherwin/webskyLowzKappa/total_zmax%.1f_zsource%.1f_kap.fits"%(zone,ztwo))
    print(hp.get_nside(mapK))
    mapFGAlm = hp.read_alm("/global/cscratch1/sd/maccrann/cmb/fg_outputs/websky_outputs/allfgs_nemo-wdr6dn_tsz-nemo-mask-snr5-mr6_ps-model-snr5/fg_terms_0/kappa_fg_" +FG+ ".fits")#kappa_fg_0093_qe.fits  
    mapFGAlm[np.isnan(mapFGAlm)]=0.0
    realMapFG = hp.alm2map(mapFGAlm,512)

    print('great')
    clsO = hp.anafast(mapK,realMapFG,lmax=1024)###cs.alm2cl(mapFG, mapK2)
    ellsO = numpy.arange(len(clsO))
    clsBinned = bin_spectrum(ellsO,clsO,ell_bin_edges)
    ellsBinned = (ell_bin_edges[1:]+ell_bin_edges[:-1])/2
    
    pickle.dump(clsBinned,open(outputDir+'biasesXcorrZ%.1f'%zone+FG+'.pkl','wb'))
    zmidarray += [zone+0.05]
    clbiasarray += [clsBinned]
    
    print('success!')
    clsRef = hp.anafast(mapK,mapCMBK,lmax=1024)###cs.alm2cl(mapFG, mapK2)
    ellsRef = numpy.arange(len(clsRef))
    clsRefBinned = bin_spectrum(ellsRef,clsRef,ell_bin_edges)
    pickle.dump(clsRefBinned,open(outputDir+'spectraXcorrZ%.1f'%zone+FG+'.pkl','wb'))
    count += 1
  else:
    zone = i*0.2+0.2
    ztwo = zone + 0.1
    zmidarray += [zone+0.05]
    clsBinned = pickle.load(open(outputDir+'biasesXcorrZ%.1f'%zone+FG+'.pkl','rb'))
    clsRefBinned = pickle.load(open(outputDir+'spectraXcorrZ%.1f'%zone+FG+'.pkl','rb'))
    print('success!')
    ellsBinned = (ell_bin_edges[1:]+ell_bin_edges[:-1])/2
    count += 1
    b[count,:] = clsBinned
    br[count,:] = clsRefBinned
  if count == 1:  
    pylab.plot(ellsBinned,clsBinned/clsRefBinned,label=str(zone),color='k')
  elif count == 10:
    pylab.plot(ellsBinned,clsBinned/clsRefBinned,label=str(zone))
  elif count == 20:
    pylab.plot(ellsBinned,clsBinned/clsRefBinned,label=str(zone))
  else:
    pylab.plot(ellsBinned,clsBinned/clsRefBinned)
 pylab.legend()
 pylab.xlim(0,1000)
 pylab.ylim(-1.5e-1,0.e-2)
 pylab.savefig(outputDir+'fractionalBiasDeltaKBins'+FG+'.png')
 pylab.clf()
 count = 0
 matrix = pickle.load(open(outputDir+'matrix.pkl','rb'))
 zvals = pickle.load(open(outputDir+'zvals.pkl','rb'))
 for i in numpy.arange(23):
    dndz = numpy.zeros(23)
    dndz[i] = 1.
    a = numpy.dot(dndz,numpy.dot(matrix,b))
    ar = numpy.dot(dndz,numpy.dot(matrix,br))
    count+=1
    if count == 1:
        pylab.plot(ellsBinned,a/ar,label=str(zvals[i]),color='k')
    elif count == 10:
        pylab.plot(ellsBinned,a/ar,label=str(zvals[i]))
    elif count == 20:
        pylab.plot(ellsBinned,a/ar,label=str(zvals[i]))
    else:
    #pylab.plot(ellsBinned,clsBinned/clsRefBinned)
        pylab.plot(ellsBinned,a/ar)
 pylab.legend(loc='best')
#pylab.ylim(0,0.05)
 pylab.savefig(outputDir+'fractionalBiasDeltaZBins'+FG+'.png')


 pylab.clf()

 Zvals = numpy.loadtxt(outputDir+'bsml_dndz_'+sample+'.txt')[:,0]#numpy.arange(100)/20.
 bdNDZ = numpy.loadtxt(outputDir+'bsml_dndz_'+sample+'.txt')[:,1]#numpy.exp(-(Zvals-1.**2.)/2./0.3**2.)*(1+Zvals)
 dndz = numpy.interp(zvals,Zvals,bdNDZ,left=0.,right=0.)
 pylab.plot(zvals,dndz)
 pylab.savefig(outputDir+sample+'dndz.png')
 pylab.clf()
 Clbias = numpy.dot(dndz,numpy.dot(matrix,b))
 Clxcorr = numpy.dot(dndz,numpy.dot(matrix,br))


 refArray = numpy.array(['ell','0093_prh', '0093_psh', '0093_qe', '0145_prh', '0145_psh', '0145_qe', 'freqcoadd_prh', 'freqcoadd_psh', 'freqcoadd_qe', 'hilc-cibd_prh', 'hilc-cibd_psh', 'hilc-cibd_qe', 'hilc-tszd_prh', 'hilc-tszd_psh', 'hilc-tszd_qe', 'hilc_prh', 'hilc_psh', 'hilc_qe'])
# print(FG)
# print(numpy.where(refArray==FG)[0])
 gerritXcorr = numpy.loadtxt(outputDir+'fg_bias_'+sample+'.dat')[:,numpy.where(refArray==FG)[0]] 
 gerritL = numpy.loadtxt(outputDir+'fg_bias_'+sample+'.dat')[:,0] 

 pylab.plot(ellsBinned,Clbias/Clxcorr,label='bias-o-matic using b dN/dz, '+FG)
 pylab.plot(gerritL,gerritXcorr,label='Gerrit unWISE '+sample+' sample reference, '+FG)
 pylab.plot(ellsBinned,Clbias*0.,color='k')

 pylab.legend()
 pylab.xlim(0,1000)
 pylab.ylim(-0.07,0.04)
 pylab.xlabel(r"$L$")
 pylab.ylabel(r"$\Delta C_L^{FG,\kappa g} / C_L^{\kappa g}$")
 pylab.savefig(outputDir+'fractionalBiasTest'+FG+''+sample+'.png',bbox_inches='tight')
 pylab.clf()
 pylab.plot(ellsBinned,Clxcorr,label='bias z=1, b=2, sigma=0.3')
 pylab.legend()
 pylab.savefig(outputDir+'xcorrSpectrumTest'+FG+''+sample+'.png')



#print(len(ellsBinned[:-1]))

