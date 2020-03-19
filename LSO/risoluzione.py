import numpy as np
import os
import logging
import argparse
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.special import erf
from scipy import signal
logging.basicConfig(level=logging.INFO)

_description = 'Analize spectrum and find energy resolution'

def change_extension(file_csv):
    """
    This function aloows you to change files' format from csv to txt.
    """
    if len(file_csv) != 0:    
        for i in range(0,len(file_grezzi),1):
            
            source = file_csv[i]
            filename, file_extension = os.path.splitext(source)
            dest = filename+'.txt'
            os.rename(source, dest)
            print('File correctly transformed from csv to txt format!')
    else:
        print('Files are already in txt format.')
    return

def gauss(x, C, mu, sigma):
    return C * norm.pdf(x, mu, sigma)

def dg(x,a1,b1,c1,a2,b2,c2):
    """
    Double gaussian distribution.
    """
    return a1*np.exp(-((x-b1)/c1)**2) + a2*np.exp(-((x-b2)/c2)**2)

def gaussian_skew(x,A,x0,sigma,alpha):    #-4,1,0
    return gauss(x,A,x0,sigma)*(1+erf(alpha*(x-x0)/sigma))

def chi2(xdata,ydata,f,*popt):
    """
    This function returns chi squared value.
    """
    mask = ydata > 0
    chi2 = sum(((ydata[mask] - f(xdata[mask], *popt)) / np.sqrt(ydata[mask]))**2.)
    nu = mask.sum() - len(popt)
    print('popt',popt)
    print('ydata e nu e popt',len(ydata),nu, len(popt))
    sigma = np.sqrt(2 * nu)
    print('chi_square     = {}'.format(chi2))
    print('expected value = {} +/- {}'.format(nu,sigma))
    return chi2

def linear(x1,x2,y1,y2,item):
    """
    Find angular coefficient and intercept of a straight line given two points.
    """
    m = (y2-y1)/(x2-x1)
    q = ((y1*x2) - (x1*y2))/(x2-x1)
    if args.show is not None:

        x_fit_covell = np.linspace(x1,x2,100)
        y_fit_covell = m*x_fit_covell+q
        plt.figure()
        plt.title('Retta di Covell. {}'.format(item))
        ydata,edges,_ = plt.hist(data_0,bins=100,label='Histo')
        plt.plot(x_fit_covell,y_fit_covell,label='retta di Covell')
        plt.xlabel('Energia [a.u]')
        plt.ylabel('Eventi [a.u]')
        plt.grid()
        plt.legend()
        plt.show()
    return m,q

def fit(xdata,ydata,first_extreme,second_extreme,parametri,_,item):
    
        
    _b=np.where(xdata>first_extreme)
    _c=np.where(xdata>second_extreme)
    xdata_fit=xdata[_b[0][0]:_c[0][0]-1]
    ydata_fit=ydata[_b[0][0]:_c[0][0]-1]
    x_fit=np.linspace(xdata_fit[0]-0.25*xdata_fit[0],xdata_fit[-1]+0.25*xdata_fit[-1],1000)

    if ( _ == 0):
        popt,pcov = curve_fit(gauss, xdata_fit, ydata_fit,p0=parametri)
        y_fit=gauss(x_fit,*popt)
        chiquadro=chi2(xdata_fit,ydata_fit,gauss,*popt)
        '''popt,pcov = curve_fit(gaussian_skew, xdata_fit, ydata_fit,p0=parametri)
        y_fit=gaussian_skew(x_fit,*popt)
        chiquadro=chi2(xdata_fit,ydata_fit,gaussian_skew,*popt)'''
    elif ( _ == 1):
        popt,pcov = curve_fit(gauss, xdata_fit, ydata_fit,p0=parametri)
        y_fit=gauss(x_fit,*popt)
        chiquadro=chi2(xdata_fit,ydata_fit,gauss,*popt)
        '''popt,pcov = curve_fit(gaussian_skew, xdata_fit, ydata_fit,p0=parametri)
        y_fit=gaussian_skew(x_fit,*popt)
        chiquadro=chi2(xdata_fit,ydata_fit,gaussian_skew,*popt)'''
    elif ( _ == 2):
        popt,pcov = curve_fit(dg, xdata_fit, ydata_fit,p0=parametri)
        y_fit=dg(x_fit,*popt)
        chiquadro=chi2(xdata_fit,ydata_fit,dg,*popt)
    elif ( _ == 3):
        popt,pcov = curve_fit(gauss, xdata_fit, ydata_fit,p0=parametri)
        y_fit=gauss(x_fit,*popt)
        chiquadro=chi2(xdata_fit,ydata_fit,gauss,*popt)
        '''popt,pcov = curve_fit(gaussian_skew, xdata_fit, ydata_fit,p0=parametri)
        y_fit=gaussian_skew(x_fit,*popt)
        chiquadro=chi2(xdata_fit,ydata_fit,gaussian_skew,*popt)'''
    else :
        popt,pcov = curve_fit(dg, xdata_fit, ydata_fit,p0=parametri)
        y_fit=dg(x_fit,*popt)
        chiquadro=chi2(xdata_fit,ydata_fit,dg,*popt)
    
    if args.show is not None:
    
        plt.figure()
        plt.title('Fit iniziale per stimare i canali di Covell. {}'.format(item))
        ydata,edges,_ = plt.hist(data_0,bins=100,label='Histo')
        plt.plot(x_fit,y_fit,label='Fit')
        plt.xlabel('Energia [a.u]')
        plt.ylabel('Eventi [a.u]')
        plt.grid()
        plt.legend()
        plt.show()
    return popt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=_description)
    parser.add_argument('-s', '--show', help='Do you want to show the images?')
    args = parser.parse_args()

    f = open('results/risoluzione.txt', 'w') 
    f.write('{} \t{} \t \n '.format('Nome','Risoluzione'))
    file_grezzi = glob.glob('*.csv')
    change_extension(file_grezzi)
    files = glob.glob('*txt') 
    for _, item in enumerate(files):
        print('----- {} -----'.format(item))
        print(_)
        """
        Load the file and obtain the histogram.
        """
        data_0 = np.loadtxt(item)
        ydata,edges,__ = plt.hist(data_0,bins=100)
        xdata = 0.5 * (edges[1:] + edges[:-1])

        """
        Find the border channel manually defining parameters.
        """
        if ( _ == 0):
            parametri=(100,3.93,0.5)         #gauss
            #parametri=(100,3.93,0.5,-4)      #skew_gauss   
            first_extreme=3.21
            second_extreme=4.65
        elif ( _ == 1):
            parametri=(100,3.35,0.5)        #gauss
            #parametri=(100,3.35,0.5,-4)     #skew_gauss
            first_extreme=2.59
            second_extreme=4.27
        elif ( _ == 2):
            parametri=(3140,2.486,0.1881,3100,2.175,0.492)
            first_extreme=1.61
            second_extreme=3.09
        elif ( _ == 3):
            parametri=(100,2.76,0.5)        #gauss
            #parametri=(100,2.76,0.5,-4)     #skew_gauss
            first_extreme=2.44
            second_extreme=3.2
        else :
            parametri=(5766,2.504,0.1836,3007,2.061,0.2861)
            first_extreme=1.61
            second_extreme=3.09
        popt=fit(xdata,ydata,first_extreme,second_extreme,parametri,_,item)

        """
        Find the exact extreme and calculate Covell line.
        Sigma's coefficients change with the fit method chosen.
        """
        if ( _ == 0):
            _d=np.where(xdata>popt[1]-2*popt[2])
            _e=np.where(xdata>popt[1]+0.5*popt[2])
        elif ( _ == 1):
            _d=np.where(xdata>popt[1]-3*popt[2])
            _e=np.where(xdata>popt[1]+2*popt[2])
        elif ( _ == 2):
            _d=np.where(xdata>popt[4]-1*popt[5])
            _e=np.where(xdata>popt[1]+2*popt[2])
        elif ( _ == 3):
            _d=np.where(xdata>popt[1]-2*popt[2])
            _e=np.where(xdata>popt[1]+1*popt[2])
        else :
            _d=np.where(xdata>popt[4]-2*popt[5])
            _e=np.where(xdata>popt[1]+2*popt[2])
        m,q=linear(xdata[_d[0][0]],xdata[_e[0][0]],ydata[_d[0][0]],ydata[_e[0][0]],item)

        """
        Subtract line's points from the first histogram. Do a gaussian (or double gaussian fit) and find the energetic resolution
        """
        ydata__diff = ydata-(m*xdata+q)
        ydata_diff = np.where(ydata__diff>0,ydata__diff,0)
        popt_new = fit(xdata,ydata_diff,xdata[_d[0][0]],xdata[_e[0][0]],parametri,_,item)
        if (_ == 2 or _ == 4):
            R1=2.35*popt_new[5]/popt_new[4]
            R2=2.35*popt_new[2]/popt_new[1]
        else:
            R1 = 2.35*popt_new[2]/popt_new[1]
            R2=0
        logging.info('R = {} , {}'.format(R1,R2))

        f.write('{} \t{} \t{} \t \n '.format(item,R1,R2))
    f.close()


   
