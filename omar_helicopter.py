#%%% Omar's helicopter code

from scipy.interpolate import UnivariateSpline

import numpy as np

import scipy.signal as sig

 

# numFREQ_NAN is the number of elements that are removed before and after the detected peak

numFREQ_NAN=5

# numFREQ_INTERP is the number of elements before and after the detected peak are used for interpolation.

numFREQ_INTERP=20

#  maximum frequency that is integrated for estimating power by PSD intergration

fMAX=100

 

# ff1 is the vector of frequencies

ind100=np.where((ff1>2)&(ff1<fMAX))[0]

ff1S=ff1[ind100]

 

# valST is 3 dimension array that stores spectrograms for multiple channels

# index 0 is channel number

# index 1 is the frequency, and

# indes 2 is time.

ve=valST[0,ind100,:].copy()

 

veST=valST[0,ind100,:].copy()

 

 

 

for resi in range(ve.shape[1]):

 

    df=ff1S[2]-ff1S[1]

    valPSD=ve[:,resi].copy()

    valPSD2=ve[:,resi].copy()

    # sig.find_peaks with int(4/df) de separation between pekas to be defiend as peaks

    peaks, dw = sig.find_peaks(valPSD, distance=int(4/df),prominence=[None,None])

    pf=dw['prominences']

    lf=dw['left_bases']

    rf=dw['right_bases']

   

    #ratio filters peaks that are strong

    ratio=valPSD[peaks]/((valPSD[lf]+valPSD[rf])/2)

   

    for jj,pi in enumerate(peaks):

        # pi has the index of the peak in the spectrogram

        # I am looking for peaks  above 10Hz the helicopter is 13 Hz +/- 2 Hz

        if ff1[pi]<10:

            continue

        #posL:posR is the section in the spectrogram around the peak that is removed

        posL=pi-numFREQ_NAN #lf[jj]

        posR=pi+numFREQ_NAN #rf[jj]

 

        if ratio[jj]>2:

            ve[posL:posR,resi]=np.nan

            valPSD2[posL:posR]=np.nan

            #sec2interp is a region around the peak tha is used for the intepolation,

            #sec2interp is larger than posL:posR

            sec2interp=valPSD2[pi-numFREQ_INTERP:pi+numFREQ_INTERP]

            # the following 2 steps are need for UnivariateSpline

            w = np.isnan(sec2interp)

            sec2interp[w]=0.

            #interpolation

            spl = UnivariateSpline(ff1S[pi-numFREQ_INTERP:pi+numFREQ_INTERP], sec2interp, w=~w,k=3) 

            # interpolated area is replaced.

            valPSD2[pi-numFREQ_NAN:pi+numFREQ_NAN]=spl(ff1S[pi-numFREQ_NAN:pi+numFREQ_NAN])

    veST[:,resi]=valPSD2

# veST is the modified spectrogram.

 