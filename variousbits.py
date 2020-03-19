from pylab import *
from scipy.optimize import curve_fit
import urllib.request, zipfile, io
dataurl = 'https://docs.google.com/uc?export=download&id=183pCEe2KtC__kM9RO-hXJqFTc3qUGEsl'
zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(dataurl).read())).extractall()
files=['BGO1.F18.csv', 'BGO2.F18.csv', 'BGO3.F18.csv',
       'LSO1.F18.csv', 'LSO2.F18.csv', 'LSO3.F18.csv',
       'NaI1.Tc99m.csv', 'NaI2.Tc99m.csv', 'NaI3.Tc99m.csv',
       'LSO1.Ba133.csv', 'LSO2.Ba133.csv']
histograms = [histogram(genfromtxt(n),bins=200) for n in files]
histograms = [(h,e,(e[:-1]+e[1:])/2) for h,e in histograms] # calcolo i centri dei bin
print('Data loaded successfully.')

#trying to plot Klein-Nishina compton behaviour
r0=2.818*10**-12
def alpha(x):
  return x/(8.2*10**-14)
def pcomp(x):
  return np.pi*r0**2/alpha(x)*((1-2*(alpha(x)+1)/alpha(x)**2)*np.log(2*alpha(x)+1)+0.5+4/alpha(x)-1/(2*(2*alpha(x)+1)**2))
xc=np.linspace(0,1,100)
plt.plot(xc,pcomp(xc))

#performing gaussian fit using only right hand side of peak
f, axes = subplots(4,3,figsize=(20,10))
axes.ravel()[-1].remove()
peaks,_= find_peaks(h,prominence=200)
cutoff=[3.75,3,4.4,2.5,3.9,3.3,1.2,1.02,0.57,2.4,2.4]
for i,n in enumerate(files):
  
    


    ax = axes[i//3,i%3]
    h, e, c = histograms[i]
    ax.step(e[:-1],h,label=n)
    ax.set_title(n)
    def gaussian_with_offset(x,A,x0,sigma,offset):
      return gaussian(x,A,x0,sigma) + offset
    hc=h
    cc=c
    hc=hc[cc>cutoff[i]]
    cc=cc[cc>cutoff[i]]

    A  = hc.max()
    x0 = float(max(e[:-1][h==hc.max()]))
    HM = e[:-1][(h>A/2)*(e[:-1]>float(max(e[:-1][h==h[peaks.max()]]))*0.8)]
    sigma = (HM[-1] - HM[0])/2.355
    popt = (A, x0, sigma, 0)
    popt, pcov = curve_fit(gaussian_with_offset, cc, hc, p0=popt,maxfev=5000)
    A, x0, sigma, offset = popt
    ax.step(e[:-1],gaussian_with_offset(e[:-1],*popt),label='fit')
    ax.set_title(f'Fit estimation: A = {A:.0f}, x0 = {x0:.2f}, sigma = {sigma:.2f}');
    ax.set_xlim((0,x0+6*sigma))

    scale_factor = 511e3/x0 # 511 keV / photopeak center
    
    print(f'Photopeak center  : {x0:.2f}')
    print(f'Energy resolution : {100*sigma*2.355/x0:.1f} %')
    print( 'Fit model         : Gaussian + offset')
f.tight_layout()

#Doing it the ways you guys did it
f, axes = subplots(4,3,figsize=(20,10))
axes.ravel()[-1].remove()
from scipy.signal import find_peaks
for i,n in enumerate(files):
  

    ax = axes[i//3,i%3]
    h, e, c = histograms[i]
    ax.step(e[:-1],h,label=n)
    ax.set_title(n)
    peaks,_= find_peaks(h,prominence=200)
    print(h[peaks.max()])
    def m(x1,y1,x2,y2):
      return (y2-y1)/(x2-x1)

    def const(x1,y1,x2,y2):
      return y1-m(x1,y1,x2,y2)*x1

    def gaussian(x,A,x0,sigma):
        return A*exp(-(x-x0)**2/(2*sigma**2))/(sigma*sqrt(2*pi))

    A  = h[peaks.max()]
    x0 = float(max(e[:-1][h==h[peaks.max()]]))
    print(float(max(e[:-1][h==h[peaks.max()]]))*0.8)
    HM = e[:-1][(h>A/2)*(e[:-1]>float(max(e[:-1][h==h[peaks.max()]]))*0.8)]
    sigma = (HM[-1] - HM[0])/2.355
    x1=e[:-1][e[:-1]>(HM[0]-sigma*2)]
    x2=e[:-1][e[:-1]>(HM[-1]+sigma*2)]
    x1=x1.min()
    print(HM[-1]+sigma*3)
    print(x2)
    x2=x2.min()
    y1=float(h[e[:-1]==x1])
    y2=float(h[e[:-1]==x2])
    popt = (A,x0,sigma)
    print(popt)
    puregauss=h-(m(x1,y1,x2,y2)*e[:-1]+const(x1,y1,x2,y2))
    popt, pcov = curve_fit(gaussian, c, puregauss, p0=popt)
    print(popt)
    ax.step(e[:-1],puregauss)
    A, x0, sigma = popt
    ax.step(e[:-1],gaussian(e[:-1],*popt),label='fit')
    ax.set_title(f'Fit estimation: A = {A:.0f}, x0 = {x0:.2f}, sigma = {sigma:.2f}');
    ax.set_xlim((0,x0+6*sigma))

    scale_factor = 511e3/x0 # 511 keV / photopeak center
    
    print(f'Photopeak center  : {x0:.2f}')
    print(f'Energy resolution : {100*sigma*2.355/x0:.1f} %')
    print( 'Fit model         : Gaussian + offset')
f.tight_layout()
