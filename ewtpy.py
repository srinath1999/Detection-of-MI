import numpy as np
#%EWT functions 
def EWT1D(f, type = "fft",N = 5, log = 0,detect = "locmax", completion = 0, reg = 'average', lengthFilter = 10,sigmaFilter = 5):
    #signal spectrum
    if(type == "fft"):
        ff = np.fft.fft(f)
        ff = abs(ff[0:int(np.ceil(ff.size/2))])#one-sided magnitude
    elif(type =="fbse"):
        ff = FBSE(f)
    #extract boundaries of Fourier Segments
    boundaries = EWT_Boundaries_Detect(ff,log,detect,N,reg,lengthFilter,sigmaFilter)
    boundaries = boundaries*np.pi/round(ff.size)
    
    if completion == 1 and len(boundaries)<N-1:
        boundaries = EWT_Boundaries_Completion(boundaries,N-1) 
    #Filtering
    #extend the signal by mirroring to deal with boundaries
    ltemp = int(np.ceil(f.size/2)) #to behave the same as matlab's round
    fMirr =  np.append(np.flip(f[0:ltemp-1],axis = 0),f)  
    fMirr = np.append(fMirr,np.flip(f[-ltemp-1:-1],axis = 0))
    if type == "fft":
        ffMirr = np.fft.fft(fMirr)
    elif type == "fbse":
        ffMirr = FBSE(fMirr)
    #build the corresponding filter bank
    mfb=EWT_Meyer_FilterBank(boundaries,ffMirr.size)    
    
    #filter the signal to extract each subband
    ewt = np.zeros(mfb.shape)
    for k in range(mfb.shape[1]):
        ewt[:,k] = np.real(iFBSE(np.conj(mfb[:,k])*ffMirr))  
    ewt = ewt[ltemp-1:-ltemp,:]
    
    return ewt,  mfb ,boundaries

def EWT_Boundaries_Detect(ff,log,detect, N, reg, lengthFilter,sigmaFilter):
    from scipy.ndimage.filters import gaussian_filter
    #apply log if needed
    if log == 1:
        ff = np.log(ff)
    
    #Global trend removal - TODO
    
    #Regularization 
    if reg == 'average':
        regFilter = np.ones(lengthFilter)/lengthFilter 
        presig = np.convolve(ff,regFilter,mode = 'same') #for even lenght, numpy's convolve is shifted when compared with MATLAB's
   
    elif reg == 'gaussian':
        regFilter = np.zeros(lengthFilter)
        regFilter[regFilter.size//2] = 1 #prefer odd filter lengths - otherwise the gaussian is skewed
        presig = np.convolve(ff,gaussian_filter(regFilter,sigmaFilter),mode = 'same') 
    else:
        presig = ff

    #Boundaries detection
    if detect  == "locmax":#Mid-point between two consecutive local maxima computed on the regularized spectrum
        boundaries = LocalMax(presig,N)
        
    elif detect == "locmaxmin":#extract the lowest local minima between two selected local maxima
        boundaries =  LocalMaxMin(presig,N)
        
    elif detect == "locmaxminf":#We extract the lowest local minima on the original spectrum between 
                                #two local maxima selected on the regularized signal
        boundaries =  LocalMaxMin(presig,N,fm = ff)
        
        
    #elif detect == "adaptivereg": #TODO
    
    return boundaries+1

        
        
def LocalMax(ff, N):
    N=N-1
    locmax = np.zeros(ff.size)
    locmin = max(ff)*np.ones(ff.size)
    for i in np.arange(1,ff.size-1):
        if ff[i-1]<ff[i] and ff[i]>ff[i+1]:
            locmax[i] = ff[i]
        if ff[i-1]> ff[i] and ff[i] <= ff[i+1]:
            locmin[i] = ff[i]
    N = min(N,locmax.size) 
    #keep the N-th highest maxima
    maxidxs = np.sort(locmax.argsort()[::-1][:N])
    #middle point between consecutive maxima
    bound = np.zeros(N)
    for i in range(N):
        if i == 0:
            a = 0
        else:
            a = maxidxs[i-1]
        bound[i] = (a + maxidxs[i])/2
        
    return bound

def LocalMaxMin(f,N,fm = 0):
    
    locmax= np.zeros(f.size)  
    if type(fm) == int: 
        f2 = f
    else:
        f2 = fm
    locmin = max(f2)*np.ones(f2.size)
    #detect local minima and maxima
    for i in np.arange(1,f.size-1):
        if ((f[i-1]<f[i]) and (f[i]>f[i+1])):
            locmax[i] = f[i]
        if ((f2[i-1]>f2[i]) and (f2[i]<f2[i+1])):
            locmin[i] = f2[i]
    #keep the N-th highest maxima and their index
    if N!=-1:
        N = N-1
        #keep the N-th highest maxima
        Imax = np.sort(locmax.argsort()[::-1][:N])
        #detect the lowest minima between two consecutive maxima
        bound = np.zeros(N)
        for i in range(N):
            if i == 0:
                a = 1
            else:
                a = Imax[i-1]
            lmin = np.sort(locmin[a:Imax[i]])
            ind = np.argsort(locmin[a:Imax[i]])
            tmpp = lmin[0]
            n = 0
            if n < len(lmin):
                n = 1
                while ((n<len(lmin)) and (tmpp==lmin[n])):
                    n = n+1
            bound[i] = a + ind[n//2]-1
    else:
        k = 0
        for i in range(locmin):
            if locmin[i]<max(f2):
                bound[k] = i-1
                k = k+1
    return bound
        
#TODO def Adaptive_Bounds_Adapt():

    
def EWT_Boundaries_Completion(boundaries,NT):
    Nd=NT-len(boundaries)
    deltaw=(np.pi-boundaries[-1])/(Nd+1)
    for k in range(Nd):
        boundaries = np.append(boundaries,boundaries[-1]+deltaw)

def EWT_Meyer_FilterBank(boundaries,Nsig):
    Npic = len(boundaries)
    #compute gamma
    gamma = 1
    for k in range(Npic-1):
        r = (boundaries[k+1]-boundaries[k])/ (boundaries[k+1]+boundaries[k])
        if r < gamma:
            gamma = r
    r = (np.pi - boundaries[Npic-1])/(np.pi + boundaries[Npic-1])
    if r <gamma:
        gamma = r
    gamma = (1-1/Nsig)*gamma#this ensure that gamma is chosen as strictly less than the min

    
    mfb = np.zeros([Nsig,Npic+1])

    #EWT_Meyer_Scaling
    Mi=int(np.floor(Nsig/2))
    w=np.fft.fftshift(np.linspace(0,2*np.pi - 2*np.pi/Nsig,num = Nsig))
    w[0:Mi]=-2*np.pi+w[0:Mi]
    aw=abs(w)
    yms=np.zeros(Nsig)
    an=1./(2*gamma*boundaries[0])
    pbn=(1.+gamma)*boundaries[0]
    mbn=(1.-gamma)*boundaries[0]
    for k in range(Nsig):
       if aw[k]<=mbn:
           yms[k]=1
       elif ((aw[k]>=mbn) and (aw[k]<=pbn)):
           yms[k]=np.cos(np.pi*EWT_beta(an*(aw[k]-mbn))/2)
    yms=np.fft.ifftshift(yms) 
    mfb[:,0] = yms
    
    #generate rest of the wavelets
    for k in range(Npic-1):
        mfb[:,k+1] = EWT_Meyer_Wavelet(boundaries[k],boundaries[k+1],gamma,Nsig)

    mfb[:,Npic] = EWT_Meyer_Wavelet(boundaries[Npic-1],np.pi,gamma,Nsig)
    
    return mfb


def EWT_beta(x):
    if x<0:
        bm=0
    elif x>1:
        bm=1
    else:
        bm=(x**4)*(35.-84.*x+70.*(x**2)-20.*(x**3))
    return bm

def EWT_Meyer_Wavelet(wn,wm,gamma,Nsig):
    Mi=int(np.floor(Nsig/2))
    w=np.fft.fftshift(np.linspace(0,2*np.pi - 2*np.pi/Nsig,num = Nsig))
    w[0:Mi]=-2*np.pi+w[0:Mi]
    aw=abs(w)
    ymw=np.zeros(Nsig)
    an=1./(2*gamma*wn)
    am=1./(2*gamma*wm)
    pbn=(1.+gamma)*wn
    mbn=(1.-gamma)*wn
    pbm=(1.+gamma)*wm
    mbm=(1.-gamma)*wm

    for k in range(Nsig):
       if ((aw[k]>=pbn) and (aw[k]<=mbm)):
           ymw[k]=1
       elif ((aw[k]>=mbm) and (aw[k]<=pbm)):
           ymw[k]=np.cos(np.pi*EWT_beta(am*(aw[k]-mbm))/2)
       elif ((aw[k]>=mbn) and (aw[k]<=pbn)):
           ymw[k]=np.sin(np.pi*EWT_beta(an*(aw[k]-mbn))/2)

    ymw=np.fft.ifftshift(ymw)
    return ymw
            

def FBSE(f):
    fbse=np.zeros(len(f))
    lin=np.linspace(np.pi/len(f),np.pi,len(f))
    for x in range(1,len(f)):#this for loop is for fbse same like dft 
        fbse = fbse + 2*lin*f[x-1]*np.sin(lin*(x-1))
    return fbse

def iFBSE(ff):
    infbse = np.zeros(len(ff))
    lin = np.linspace(0.0,np.pi-np.pi/len(ff),len(ff))
    for x in range(1,len(ff)):
        infbse = infbse + (ff[x-1]/(x*np.pi))*np.sin(lin*x)
    return infbse





