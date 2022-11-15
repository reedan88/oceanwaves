# +
import gsw
import numpy as np
from scipy import signal, fft
import matplotlib.pyplot as plt


def wavenumber(f, h, lat):
    """
    Calculate the wavenumber from the wave frequency and water depth
    
    This utilizes the Hunt (1979) approximation as implemented by 
    George Voulgaris (1992)
    
    Parameters
    ----------
    f: int
        Wave frequency = 1/T (T is wave period)
    h: float
        Water depth (m)
    lat: float
        Latitude of the deployed pressure sensor (decimal degrees)
        
    Returns
    -------
    k: float
        Wave number k calculated using a polynomial approximation
    """
    # Location specific gravity
    g = gsw.grav(lat, h)
    # Radian frequency [rad/sec]
    w = 2*np.pi*f
    # Polynomial expansion parameters
    dum1 = (w**2)*h/g
    dum2 = dum1 + (1.0 + 0.6522*dum1 + 0.4622*dum1**2 + 0.0864*dum1**4 + 0.0675*dum1**5)**(-1)
    dum3 = np.sqrt(g*h*dum2**(-1)) / f
    # Calculate the wave number
    k = 2*np.pi*dum3**(-1)

    return k


def detrendHeight(pt):
    """Remove trend from a time series
    
    Fits a straight line to a vector of values using lm(), and uses the
    regression coefficients to subtract off the linear trend from the values.
    
    Typically this is used to remove a tidal trend from a ocean 
    surface height time series before attempting to calculate statistics for
    waves, which are variations above and below the mean surface height. 
    
    Returns a series of residuals around the linear trend, in the original units
    of measurement (typically meters). 
    
    Code adapted from cran package oceanwaves (Author: Luke Miller)

    Parameters
    ----------
    pt: array_like
        A vector of numeric values to be detrended
        
    Returns
    -------
    pt: numpy.array
        A vector of the input numeric values which have been detrended
    h: float
        The mean height of the seasurface
    n: int
        Length of the dataset
    (m, c): tuple(float, float)
        The slope and intercept of the linear regression        
    """
    
    # Determine length of pt
    n = len(pt)
    
    # Make a sequence of indices
    x = np.arange(0, n, 1)
    
    # Fit a simple linear regression through the segment of sea surface
    # height data and extract the intercept + slope
    x = np.arange(0, n, 1)
    A = np.vstack([x, np.ones(n)]).T
    m, c = np.linalg.lstsq(A, pt, rcond=None)[0]
    
    # Calculate a depth h at the midpoint of the segment of 
    # data using the regression coefficients stored in 'trend' 
    # (i.e. intercept and slope) 
    h = c + m*n/2 

    # Remove the linear trend from the segment of sea surface heights
    pt = pt - (c + (m * x)) 
    
    # Return the detrended vector, mean height, segment length, and
    # regression coefficients
    return pt, h, int(n), (m, c)


def pressureCorrection(wave_data, Fs, z, M=512, freqLim=(0.05, 0.33), plot=False):
    """
    Correct for depth attenuation of a water surface elevation pressure signal.
    
    This functions corrects the water surface elvation time series for the 
    effects of attenuation with depth.
    
    Parameters
    ----------
    wave_data: array_like
        An array of positive wave surface elevations in meters
    Fs: float or int
        Sampling frequency in Hertz
    z: float or int
        Depth of the deployed pressure sensor (positive)
    M: int, Default=512
        Length of time series segments that will be used in the detrending and
        attenuation correction operations
    freqLim: tuple, Default=(0.05, 0.33)
        Frequency range over which to resolve wave signature. Maximum frequency
        will be overwritten in the pressure correction falls below 0.0025/Fs,
        the minimum resolvable sensitivity of the SBE26 instrument
        
    Returns
    -------
    H: numpy.array
        An array of sea surface elevations [m] corrected for pressure attenuation of
        waves that occurs with depth
    """
    # Input parameters
    pt = wave_data

    # Prep the output data
    H_with_NaN = pt           # Save initial sea surface values
    not_NaN = ~np.isnan(pt)   # Find the NaNs
    pt = pt[not_NaN]          # Remove the NaNs
    pt_size = np.size(pt)     # Size of the cleaned up sea surface values
    m = len(pt)               # Length of the cleaned up surface values
    N_overlap = M/2           # Size of segment overlengths
    N = int(np.ceil(m/M)*M)   # Length of data zero-padded to nearest multiple of M
    H = np.zeros(N)           # Pre-allocate output vector
    overlap_window = signal.windows.hann(M) # Hanning window

    # Set the frequencies
    f = np.concatenate(([np.nan], np.arange(1, M/2+1, 1)*Fs/M)) # Vector of frequencies with NaN in position [0]
    min_frequency = freqLim[0]      # minimum frequency (wave period of 20 seconds)
    max_frequency = freqLim[1]      # maximum frequency (wave period of 3 seconds)

    # Set the segments
    # Step through the segments of the sea surface elevation data
    for q in np.arange(0, N-N_overlap+1, N_overlap):
        q = int(q)
        o = int(np.min([q+M, m]))
        print(str(q) + " :: " + str(o))

        # Get an individual segment data
        pt_seg = pt[q:o]
        n_seg = len(pt_seg)

        # Detrend the data
        pt_seg_detrended, h, n, (slope, intercept) = detrendHeight(pt_seg)

        # Calculate wavenumber
        K = wavenumber(f, h, lat)

        # Calculate correction factor for the spectrum
        Kpt = np.cosh(K*(h - z)) / np.cosh(K*h)

        # Find where the correction factor drops below the maximum correction allowed
        Kpt_max = 0.0025/Fs    # Maximum correction allowed is 0.0025/sampling rate
        Kpt_max_ind = np.max(np.where(Kpt >= Kpt_max))        # Get the index of the max correction factor
        Kpt_max_freq = f[Kpt_max_ind]                     # Frequency of the maximum pressure correction
        if Kpt_max_freq < max_frequency:
            max_frequency = Kpt_max_freq

        # Replace the correction values below the min frequency and above the max frequency with 1 (i.e. no correction)
        Kpt[(f < min_frequency) | (f > max_frequency)] = 1

        # Get the index of where the max frequency is
        f_max_ind = np.max(np.where(f <= max_frequency))

        # Now linearly decrease correction factor for frequencies above the max_frequency
        Kpt_lin_ind = np.arange(f_max_ind, int(np.min((f_max_ind+np.fix(len(K)/10), len(K)))), 1)
        Kpt_lin = np.arange(len(Kpt_lin_ind), 1, -1) * (Kpt[f_max_ind]-1)/len(Kpt_lin_ind) + 1
        Kpt[Kpt_lin_ind[1:]] = Kpt_lin

        # Make the second half of the correction series symmetric
        Kpt_reversed = np.flip(Kpt[np.arange(1, int(M/2), 1)])

        # Append the reversed series to the original to get a symmetric series
        Kpt = np.concatenate([Kpt, Kpt_reversed])

        # Replace the first NaN entry with 1
        Kpt[0] = 1

        # Zero pad the segment length if it is less than the overlap window
        if n_seg < M:
            pt_seg_detrended = np.concatenate((pt_seg_detrended, np.zeros(M-n_seg)))

        # Calculate the power spectrum
        P = fft.fft(pt_seg_detrended)
        # Apply correction factor
        Pcor = P / Kpt
        # Corrected sea-surface-heights
        H_seg = np.real(fft.ifft(Pcor))
        # Truncate back to original segment length
        H_seg = H_seg[0:n_seg]
        # Add the trend back in to the reconstructed series
        H_seg = H_seg + intercept + slope*np.arange(1, n+1, 1)

        # Generate the corrected sea surface elevation 
        H_windowed = H_seg*overlap_window[0:n_seg]
        H[q:o] = H[q:o] + H_seg*overlap_window[0:n_seg]
        if q == 0:
            H[0:int(np.min((N_overlap, n_seg)))] = H_seg[0:int(np.min((N_overlap, n_seg)))]
        if q+M > N and n_seg > N_overlap:
            H[q+int(N_overlap):o] = H_seg[int(N_overlap):len(H_seg)]
            
        if plot == True:
            # Plot the segment results
            fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(15, 12))

            ax[0].plot(Kpt, marker=".", color="tab:blue")
            ax[0].set_title("Pressure Correction Value")
            ax[0].grid()

            ax[1].plot(pt[q:o], marker=".", color="tab:blue", label="Original Sea Surface Height")
            ax[1].plot(H_seg, marker=".", color="tab:red", label="Corrected Sea Surface Height")
            ax[1].set_title("Segment Height")
            ax[1].grid()
            ax[1].legend()

            ax[2].plot(pt, marker=".", color="tab:blue", label="Original Sea Surface Height")
            ax[2].plot(np.arange(q, o, 1), H_windowed, marker=".", color="tab:red", label="Windowed Corrected SSH")
            ax[2].plot(H, marker=".", color="tab:green", label="Final SSH")
            ax[2].set_title("Windowed data")
            ax[2].grid()
            ax[2].legend()

    # Truncate H back into original length m
    return H[0:m]


def waveSpectra(data, Fs):
    """Calculate wave statistics from the SSH data, corrected for pressure attentuation.
    
    This function calculates the significat wave height, peak period, the average
    period using two different methods, and the estimates of spectral width for
    sea surface height data measured by a subsurface pressure sensor. The SSH data
    should be first corrected for the effect of depth attenuation of the pressure
    signal using the pressureCorrection function. Statistics are derived using the 
    Welch PSD method.
    
    Parameters
    ----------
    data: array_like
        An array of the SSH data that has been corrected for the pressure attenuation
        effect.
    Fs: float, int
        The sampling frequency
        
    Returns
    -------
    Hm0: float
        The signicant wave height [m], calculated from the spectral moment 0
    Tp: float
        The peak period [sec], calculate from teh frequency at the spectral max
    T_0_1: float
        The average period [sec] calculated using the zeroth and first spectral 
        moments following NDBC procedure
    T_0_2: float
        The average period [sec] calculated using the zeroth and second spectral
        moments
    ESP2: float
        Spectral width parameter
    ESP4: float
        Spectral width paramter
    """
    
    # -------------------------------------------------------
    # Prepare data for spectral analysis
    n_seg = 4                 # Number of segments
    m = len(data)             # Length of the data set
    M = np.fix(m/n_seg/2)*2   # Length of the segment
        
    # Detrend the SSH data
    detrendedSSH, h, n, (slope, intercept) = detrendHeight(data)
    
    # -------------------------------------------------------
    # Calculate the spectrogram
    f, P = signal.welch(detrendedSSH, Fs, nperseg=M, noverlap=M/2)

    # Remove the zero-frequency entries
    f = f[1:]
    P = P[1:]

    # Normalize the PSD
    p = P*2/Fs

    # --------------------------------------------------------
    # Calculation of wave parameters
    integmin = np.min(np.where(f >= 0))
    integmax = np.max(np.where(f <= max_frequency * 1.5))

    # Bandwidth
    deltaf = f[1] - f[0]

    # Calculate the moments of the spectrum, from -2nd to 0th to 4th
    # For a spectrum, the 0th moment represents the variance of the data
    moment = np.zeros(7)

    for i in np.arange(-2, 5, 1):
        moment[i+2] = np.sum( (f[integmin:integmax+1]**i) * (p[integmin:integmax+1])) * deltaf

    # Peak period, calculated from the frequency at maximum of spectrum
    Tp = f[p == np.max(p)]

    # Estimated variance of the time series
    m0 = moment[2]

    # Estimate significant wave height based on the spectral moment 0
    # This value is approximately equal to the average of the highest
    # 1/3 of the waves in the time series
    Hm0 = 4 * np.sqrt(m0)

    # Calculate average period m0/m1, units seconds, following NDBC method for average period
    T_0_1 = moment[2] / moment[1+2]

    # Calculate average period (m0/m2)**(1/2), units seconds, following SIO method
    T_0_2 = (moment[2] / moment[2+2])**(0.5)

    # Spectral width estimates
    ESP2 = (moment[2] * moment[4] / moment[3]**2 - 1)**(0.5)
    ESP4 = (1 - moment[4]**2 / (moment[2] * moment[6]))**(0.5)
   
    return Hm0, float(1/Tp), T_0_1, T_0_2, ESP2, ESP4
