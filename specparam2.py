# specparam with flattened array epr channel saved 
import os
from glob import glob
import pandas as pd
from scipy.io import loadmat, savemat
from specparam import SpectralModel
import numpy as np

# set paths 
inpath = r"D:\Linus\rsEEG\Analysis2\PSD_han2s50"
mat_files = glob(os.path.join(inpath, '*.mat'))                     # remove [:1] to process all 

outpath_specparam = r"D:\Linus\rsEEG\Analysis2\SpecParam_2"
os.makedirs(outpath_specparam, exist_ok=True)

###### to change fitting range: 
upper_end = 45 
lower_end = 1

# Define frequency bands (Hz)
BANDS = {
    'Delta': (lower_end, 4),
    'Theta': (4, 8),
    'Alpha': (8, 13),
    'Beta':  (13, 30),
    'Gamma': (30, upper_end)
}

# Helper function to find the biggest peak in a specific frequency band
def get_biggest_peak_in_band(peaks, f_range):
    """Returns [CF, PW, BW] of the highest-power peak in the range, or NaNs if none exist."""
    if len(peaks) == 0:
        return np.array([np.nan, np.nan, np.nan])
    
    # Filter peaks where Center Frequency (column 0) is within f_range
    band_peaks = peaks[(peaks[:, 0] >= f_range[0]) & (peaks[:, 0] <= f_range[1])]
    
    if len(band_peaks) == 0:
        return np.array([np.nan, np.nan, np.nan])
        
    # Find the index of the peak with the maximum Power (column 1)
    max_pw_idx = np.argmax(band_peaks[:, 1])
    return band_peaks[max_pw_idx, :]

# Initialize a list to hold the global results for the CSV
global_results_list = []

for file in mat_files:
    print(f"Processing file: {file}")
    data = loadmat(file)
    chanlabels = data['chanlabels'].squeeze()
    total_windows = data['total_windows'].squeeze()
    freqs = data['freqs'].squeeze()
    psd = data['psd_avg'].squeeze()

    # channelwise SpecParam
    fms = []
    flattened_ch = []
    
    # Pre-calculate the frequency mask to align with specparam's internal fitting range
    freq_mask = (freqs >= lower_end) & (freqs <= upper_end)
    
    for ch_idx, chan in enumerate(chanlabels):
        fm = SpectralModel(peak_width_limits=[1, 18], 
                           max_n_peaks=6, 
                           min_peak_height=0.05, 
                           peak_threshold=2.0, 
                           aperiodic_mode='fixed')
        
        fm.fit(freqs, psd[ch_idx, :], [lower_end, upper_end])
        fms.append(fm)
        
        # Calculate flattened array safely: log10(power) - aperiodic_fit
        power_log = np.log10(psd[ch_idx, freq_mask])
        ap_fit = fm.get_data('aperiodic')
        flattened_ch.append(power_log - ap_fit)

    # Convert to 2D numpy array (Channels x Frequencies) for clean MATLAB export
    flattened_ch = np.array(flattened_ch) 

    # Global Specparam 
    psd_global = psd.mean(axis=0)  # average across channels
    fm_global = SpectralModel(peak_width_limits=[1, 18], 
                              max_n_peaks=6,
                              min_peak_height=0.05,
                              peak_threshold=2.0,
                              aperiodic_mode='fixed')  

    fm_global.fit(freqs, psd_global, [lower_end, upper_end])

    # Calculate global flattened spectrum for consistency
    power_log_global = np.log10(psd_global[freq_mask])
    ap_fit_global = fm_global.get_data('aperiodic')
    flattened_global = power_log_global - ap_fit_global

    # extract global peaks per band
    peaks = fm_global.get_params('peak')
    
    peak_delta = get_biggest_peak_in_band(peaks, BANDS['Delta'])
    peak_theta = get_biggest_peak_in_band(peaks, BANDS['Theta'])
    peak_alpha = get_biggest_peak_in_band(peaks, BANDS['Alpha'])
    peak_beta = get_biggest_peak_in_band(peaks, BANDS['Beta'])
    peak_gamma = get_biggest_peak_in_band(peaks, BANDS['Gamma'])

    # data extraction for CSV
    subject_id = os.path.basename(file).replace('_PSD.mat', '')
    
    global_results_list.append({
        'Filename': subject_id,
        'Offset': fm_global.get_params('aperiodic', 'offset'),
        'Exponent': fm_global.get_params('aperiodic', 'exponent'),
        'R_Squared': fm_global.get_metrics('gof_rsquared'),
        'Error_MAE': fm_global.get_metrics('error_mae'),
        'N_Peaks': len(peaks),
        # Band-specific peaks (Center Freq, Power, Bandwidth)
        'Delta_CF': peak_delta[0], 'Delta_PW': peak_delta[1], 'Delta_BW': peak_delta[2],
        'Theta_CF': peak_theta[0], 'Theta_PW': peak_theta[1], 'Theta_BW': peak_theta[2],
        'Alpha_CF': peak_alpha[0], 'Alpha_PW': peak_alpha[1], 'Alpha_BW': peak_alpha[2],
        'Beta_CF':  peak_beta[0],  'Beta_PW':  peak_beta[1],  'Beta_BW':  peak_beta[2],
        'Gamma_CF': peak_gamma[0], 'Gamma_PW': peak_gamma[1], 'Gamma_BW': peak_gamma[2],
        # Full array dump just in case
        'All_Peaks_CF_PW_BW': str(peaks.tolist()) if len(peaks) > 0 else "[]" 
    })

    # saving indiv .mat files
    out_name = os.path.basename(file).replace('_PSD.mat' , '_specparam.mat')
    out_filepath = os.path.join(outpath_specparam, out_name)
    
    savemat(out_filepath, {
        # Global parameters       
        'offset': fm_global.get_params('aperiodic', 'offset'),
        'exponent': fm_global.get_params('aperiodic', 'exponent'),
        'peaks': peaks,
        'r_squared': fm_global.get_metrics('gof_rsquared'),
        'error': fm_global.get_metrics('error_mae'),
        
        # Band-specific biggest peaks [CF, PW, BW]
        'global_peak_delta': peak_delta,
        'global_peak_theta': peak_theta,
        'global_peak_alpha': peak_alpha,
        'global_peak_beta': peak_beta,
        'global_peak_gamma': peak_gamma,

        # For plotting
        'power_spectrum': psd_global,
        'frequencies': freqs,
        'frequencies_fitted': freqs[freq_mask], # Necessary x-axis for plotting fitted/flattened data
        'aperiodic_fit': fm_global.get_data('aperiodic'),
        'periodic_fit': fm_global.get_data('peak'),
        'flattened_global': flattened_global,

        # Channelwise parameters
        'chanlabels': chanlabels,
        'flattened_ch': flattened_ch,
        'offset_ch': [fm.get_params('aperiodic', 'offset') for fm in fms],
        'exponent_ch': [fm.get_params('aperiodic', 'exponent') for fm in fms],
        'peaks_ch': [fm.get_params('peak') for fm in fms],
        'r_squared_ch': [fm.get_metrics('gof_rsquared') for fm in fms],
        'error_ch': [fm.get_metrics('error_mae') for fm in fms]
    })

    print(f"Successfully saved: {out_name}")

# save CSV
df_global = pd.DataFrame(global_results_list)
csv_outpath = os.path.join(outpath_specparam, "global_specparam_results.csv")
df_global.to_csv(csv_outpath, index=False)

print(f"\n--- ALL DONE ---")
print(f"Global results successfully saved to: {csv_outpath}")