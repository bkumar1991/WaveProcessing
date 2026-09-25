import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from pathlib import Path
import os 
import shutil
import glob
from obspy import read, Stream
from scipy.signal import hilbert
import numpy as np
from scipy.optimize import minimize
import obspy
root_path = Path('D:/python codes/Example_Data')
new_root_path = Path('D:/python codes/Example_Data_Processed')            # output path to save the final processed waveforms


def calculate_envelope_peak(trace):
    
    analytic_signal = hilbert(trace.data)
    envelope = np.abs(analytic_signal)  # Envelope is the absolute value of the analytic signal
    peak_envelope = np.max(envelope)  
    
        
    return peak_envelope
    
def calculate_envelope_peak_time(trace):
    
    analytic_signal = hilbert(trace.data)
    envelope = np.abs(analytic_signal)  # Envelope is the absolute value of the analytic signal
    peak_time = np.argmax(envelope) * trace.stats.delta  # Time of the peak (index of max envelope * delta time)
    
    return peak_time

def perform_matrix_inversion(trace_offsets, peak_envelopes, degree=1):
    
    
    G = np.vstack([trace_offsets**i for i in range(degree+1)]).T  
    
    # Perform the matrix inversion: m = (G^T * G)^-1 * G^T * y
    G_T = G.T
    try:
        
        m = np.dot(np.dot(inv(np.dot(G_T, G)), G_T), peak_envelopes)
        
        
    except np.linalg.LinAlgError as e:
        print("Error in matrix inversion: ", e)
        return None, None
    
    return m, G

def perform_matrix_inversion_for_peak_time(trace_offsets, peak_time, degree=1):
    
    # Construct the design matrix G for polynomial fitting
    G1 = np.vstack([trace_offsets**i for i in range(degree+1)]).T  # G = [1, x, x^2, ..., x^degree]
    
    # Find the indices of maximum and minimum peak times
    max_peak_indices = np.argsort(peak_time)[-2:]  # Select the indices of maximum peak times
    min_peak_indices = np.argsort(peak_time)[:2]  # Select the indices of minimum peak times
    
    # Weights: initially, all weights are 1, but we will increase the weights for max and min peak time points
    weights = np.ones(len(peak_time))
    
    weights[max_peak_indices] = 10  
     
    
    # Create a diagonal weight matrix (for weighted least squares)
    W = np.diag(weights)
    
    # Perform the matrix inversion with weighted least squares
    G1_T_W = np.dot(G1.T, W)
    
    try:
        
        m1 = np.dot(np.dot(inv(np.dot(G1_T_W, G1)), G1_T_W), peak_time)
        
    except np.linalg.LinAlgError as e:
        print("Error in matrix inversion: ", e)
        return None, None
    
    return m1, G1


def model_function(trace_offsets, m, degree=1):
    
    model_values = np.sum([m[i] * trace_offsets**i for i in range(degree+1)], axis=0)
        
    return model_values

def model_function_predicted_peak_time(trace_offsets, m1, degree=1):
    
    model_time_values = np.sum([m1[i] * trace_offsets**i for i in range(degree+1)], axis=0)
    
    return model_time_values

def predict_envelope_from_model(trace_offsets, m, degree=1):
    
    predicted_envelope = model_function(trace_offsets, m, degree)
    
    return predicted_envelope

def predict_envelope_peak_time_from_model(trace_offsets, m1, degree=1):
    
    predicted_peak_time = model_function_predicted_peak_time(trace_offsets, m1, degree)
    
    return predicted_peak_time


def save_plots_with_optimized_inversion(folder, save_path1, save_path2, save_path3, save_path4, save_path6, save_path7, save_path8, save_path9,save_path10 ):
    sac_files = glob.glob(f"{folder}/*.SAC")
        

    if not sac_files:
        return

    trace_offsets = []
    peak_envelopes = []  # List to store peak envelope values
    all_sac_files = []  # List to store the corresponding SAC file paths
    predicted_peak_time_penalty= []
    peak_times= []
    predicted_peak_times= []
    corrected_peak_times= []
    difference_peak_time= []
    polarity_reversed= []
    normal=[]
    outlier_polarity= []
    traces=[]
    correlations =[]
    correlations = []
    normal_offsets = []
    normal_corrs = []
    rt= Stream()
    nt= Stream()
    reversed_offsets = []
    reversed_corrs = []
    all_st_4 = Stream()
    
    for i, sac_file in enumerate(sac_files):
        file_path = os.path.join(folder, sac_file)
        st = read(file_path)
        trace = st[0]
        
        trace_data = trace.data
        traces.append(trace_data)

               
        trace.stats.distance = trace.stats.sac.gcarc * 110 *1000   # offset distances if available
        trace_offset = trace.stats.distance /1e3  

        
        
        
        # Calculate peak value of envelope (observed envelope)
        peak_envelope = calculate_envelope_peak(trace)
        peak_envelopes.append(peak_envelope)

        peak_time = calculate_envelope_peak_time(trace)
        peak_times.append(peak_time)
        
                
        analytic_signal = hilbert(trace.data)  # Compute the analytic signal using Hilbert transform
        envelope = np.abs(analytic_signal)  # Envelope is the absolute value of the analytic signal
        

        trace_offsets.append(trace_offset)
        all_sac_files.append(sac_file)  # Keep track of the SAC files

        

    trace_offsets = np.array(trace_offsets)
    traces.append(trace.data)
    max_len = max(len(t) for t in traces)
    #traces = [t[:max_len] for t in traces]
    traces_np = np.array(t[:max_len] for t in traces)

       
    
    m, G = perform_matrix_inversion(trace_offsets, peak_envelopes, degree=1)
    
    if m is None:
        print("Error in performing matrix inversion. Exiting function.")
        return
    m1, G1 = perform_matrix_inversion_for_peak_time(trace_offsets, peak_times, degree=1)
    
    if m1 is None:
        print("Error in performing matrix inversion. Exiting function.")
        return
    

      

    predicted_envelope_peak_time = predict_envelope_peak_time_from_model(trace_offsets, m1, degree=1)
    peak_times=np.array(peak_times)
    predicted_envelope_peak_time= np.array(predicted_envelope_peak_time)
    

    residuals = (np.array(peak_times) - predicted_envelope_peak_time)
    # print(residuals)
    residual_mean = np.mean(residuals)
    print("residual mean is :", np.abs(residual_mean))
    residual_std = np.std(residuals)
    print("residual standard deviation is :", residual_std)

    # residuals = np.array(peak_times) - predicted_envelope_peak_time

    # n = len(residuals)
    # p = 2   # number of fitted parameters

    # residual_se = np.sqrt(np.sum((residuals - np.mean(residuals))**2) / (n - p))

    # ratio = residual_se / np.std(peak_times)
    # print(ratio)

    mad = np.median(np.abs(residuals - np.median(residuals))) + 1e-10
    robust_std = 1.4826 * mad  

    if residual_std < 1.5 * robust_std:
        outlier_threshold = 3.0
    else:
        outlier_threshold = 2.3

    
    residual_outlier= np.abs(residuals - residual_mean) > outlier_threshold * residual_std
    print("residual outlier std", np.std(residual_outlier))
    residual_outlier = residual_outlier.astype(bool)
    outlier_mask = residual_outlier
    
    outlier_sac_files = [sac_file for i, sac_file in enumerate(all_sac_files) if outlier_mask[i]]
    trace_offsets_filtered = trace_offsets[~outlier_mask]
    peak_times_filtered = peak_times[~outlier_mask]
    predicted_envelope_peak_time_filtered = predicted_envelope_peak_time[~outlier_mask]

    
       
    plt.figure(figsize=(4.0, 2.3))
    plt.plot(trace_offsets_filtered, peak_times_filtered, 'ko',markersize=3.2, label= 'Observed Peak Times ')
    plt.gca().invert_yaxis()
    plt.plot(trace_offsets_filtered, predicted_envelope_peak_time_filtered, 'k-',markersize=3.2, linewidth= 1.0, label= 'Predicted Peak Times ')
    plt.gca().invert_yaxis()
    plt.scatter(trace_offsets[outlier_mask], peak_times[outlier_mask], color='red',s=6, label='Outliers', zorder=5)
    plt.gca().invert_yaxis()
    plt.xlabel('Epicentral Distance (Km)', fontsize=10)
    plt.ylabel('Envelope Peak Time (s)', fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=9)
    # plt.title('Observed vs Predicted Envelope Peak Times', fontsize=10)
    # plt.legend(fontsize=8, loc="upper right", frameon= False)
    # plt.text(0.02, 1.065, '(a)',transform=plt.gca().transAxes, fontsize=10, va='top', ha='left')
    plt.tick_params(axis='both', which='major', labelsize=5)
    plt.tight_layout()
    plt.savefig(save_path6, format="pdf",dpi=300, bbox_inches='tight' )
    plt.close()

    
    all_st_1 = Stream()
    fig3, ax3= plt.subplots(figsize=(6.5, 4))
    ax3 = fig3.axes[0]
    for i, sac_file in enumerate(sac_files):
        file_path = os.path.join(folder, sac_file)
        st = read(file_path)
        trace = st[0]
        time = np.arange(0, trace.stats.npts / trace.stats.sampling_rate, trace.stats.delta)
        trace.stats.distance = trace.stats.sac.gcarc * 110  
        # trace.stats.distance = trace.stats.sac.gcarc * 110 *1000   # offset distances if available
        # trace_offset = trace.stats.distance /1e3  

        
        all_st_1 += trace 
        
    all_st_1.plot(type='section', fig=fig3, time_down=True, scale=2.5, outfile=save_path7)
    ax3.set_ylim(1000,200)
    ax3.set_ylabel('Time (s)', fontsize=10)
    ax3.set_xlabel("Epicentral Distance (Km)", fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=9)
    # plt.text(0.02, 1.04, '(b)',transform=plt.gca().transAxes, fontsize=10, va='top', ha='left')
    plt.tight_layout()
    plt.savefig(save_path7)
    
    # Delete outlier files
    
    for i, sac_file in enumerate(outlier_sac_files):
        outlier_file_path = os.path.join(folder, sac_file)
                
        if os.path.exists(outlier_file_path):
            print(f"Deleting residual file: {outlier_file_path}")
            os.remove(outlier_file_path)

    sac_files = glob.glob(f"{folder}/*.SAC")
    fig2 = plt.figure(figsize=(6.5, 4))
    all_st_3 = Stream()
    for i, sac_file in enumerate(sac_files):

        file_path = os.path.join(folder, sac_file)
        st = read(file_path)
        trace = st[0]
        
        time = np.arange(0, trace.stats.npts / trace.stats.sampling_rate, trace.stats.delta)
        trace.stats.distance = trace.stats.sac.gcarc * 110 
        # trace.stats.distance = trace.stats.sac.gcarc * 110 *1000 
        # trace_offset = trace.stats.distance /1e3
        all_st_3 += trace 
        
    all_st_3.plot(type='section', fig=fig2, time_down=True, scale=2.5, outfile=save_path4)
    ax2 = fig2.axes[0]
    ax2.set_ylim(1000, 200)
    ax2.set_xlabel("Epicentral Distance (Km)", fontsize=10)
    ax2.set_ylabel('Time (s)', fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=9)
    plt.tight_layout()
    plt.savefig(save_path9)
    plt.close(fig2)
    
    sac_files = glob.glob(f"{folder}/*.SAC")

    # Reset everything
    traces = []
    trace_offsets = []
    peak_times = []
    peak_envelopes = []
    all_sac_files = []

    for sac_file in sac_files:

        file_path = sac_file   
        st = read(file_path)
        trace = st[0]

        trace_data = trace.data
        traces.append(trace_data)

        trace.stats.distance = trace.stats.sac.gcarc * 110 * 1000
        trace_offset = trace.stats.distance / 1e3

        peak_envelope = calculate_envelope_peak(trace)
        peak_time = calculate_envelope_peak_time(trace)

        trace_offsets.append(trace_offset)
        peak_envelopes.append(peak_envelope)
        peak_times.append(peak_time)

        all_sac_files.append(os.path.basename(sac_file))

    # Convert to numpy
    trace_offsets = np.array(trace_offsets)
    peak_times = np.array(peak_times)
    traces.append(trace_data)

         
    max_len = max(len(t) for t in traces)
    # traces = [t[:max_len] for t in traces]
    traces_np = np.array([t[:max_len] for t in traces])
    min_len = min(len(traces_np), len(all_sac_files), len(trace_offsets))
       
       
    with open(output_txt_file, 'w') as f1:
                
        
        reference_trace = None
        reference_trace_num = None
        prev_reference_num = -1
        coherence_threshold = 0.5
        high_corr_threshold = 0.7
        neg_streak_threshold = 2   
        
        neg_streak = 0

        for i in range(min_len):
        
            current_trace = traces_np[i]

              
            if i == 0:                       # neighbor correlations
                prev_corr = 0
            else:
                prev_corr = np.corrcoef(current_trace, traces_np[i-1])[0, 1]

            if i == min_len - 1:
                next_corr = 0
            else:
                next_corr = np.corrcoef(current_trace, traces_np[i+1])[0, 1]

            in_trend = (prev_corr >= coherence_threshold or next_corr >= coherence_threshold)

            if reference_trace is None:

                reference_trace = current_trace
                reference_trace_num = i + 1

                neg_streak = 0

                normal.append(all_sac_files[i])
                normal_offsets.append(trace_offsets[i])
                normal_corrs.append(1.0)

                f1.write(
                    f"{'trace':<6}: {i+1:<5}  "
                    f"{'reference trace':<16}: {reference_trace_num:<5} "
                    f"{'corr_normal':<12}: {1.0:<8.4f}  "
                    f"{'offset':<10}: {trace_offsets[i]:<12.4f}\n"
                )
                continue
            
            if not in_trend:
                corr_ref = np.corrcoef(current_trace, reference_trace)[0, 1]

                if corr_ref < 0:
                    neg_streak +=1    
                    offset_jump = ( i >0 and
                        abs(trace_offsets[i] - trace_offsets[i-1]) > 50)

                    if offset_jump:

                        prev_reference_num = reference_trace_num
                        reference_trace = current_trace.copy()
                        reference_trace_num = i + 1

                        neg_streak = 0        
                        normal.append(all_sac_files[i])
                        normal_offsets.append(trace_offsets[i])
                        normal_corrs.append(abs(corr_ref))

                        f1.write(
                            f"{'trace':<6}: {i+1:<5}  "
                            f"{'reference trace':<16}: {prev_reference_num:<5}  "
                            f"{'corr_normal':<12}: {abs(corr_ref):<8.4f}  "
                            f"{'offset':<10}: {trace_offsets[i]:<12.4f}  "
                            
                        )
                    elif neg_streak >= 2:

                        normal.append(all_sac_files[i])
                        normal_offsets.append(trace_offsets[i])
                        normal_corrs.append(abs(corr_ref))

                        f1.write(
                            f"{'trace':<6}: {i+1:<5}  "
                            f"{'reference trace':<16}: {prev_reference_num:<5}  "
                            f"{'corr_normal':<12}: {abs(corr_ref):<8.4f}  "
                            f"{'offset':<10}: {trace_offsets[i]:<12.4f}  "
                            
                        )     

                    else:

                        polarity_reversed.append(all_sac_files[i])
                        reversed_offsets.append(trace_offsets[i])
                        reversed_corrs.append(corr_ref)

                        f1.write(
                            f"{'trace':<6}: {i+1:<5}  "
                            f"{'reference trace':<16}: {reference_trace_num:<5} "
                            f"{'corr_normal':<12}: {corr_ref:<8.4f}  "
                            f"{'offset':<10}: {trace_offsets[i]:<12.4f}"     
                            f"REMOVED\n"
                        )
                
                                
                else:
                    prev_reference_num = reference_trace_num
                    normal.append(all_sac_files[i])
                    normal_offsets.append(trace_offsets[i])
                    normal_corrs.append(corr_ref)

                    f1.write(
                        f"{'trace':<6}: {i+1:<5}  "
                        f"{'reference trace':<16}: {prev_reference_num:<5} "
                        f"{'corr_normal':<12}: {1.0:<8.4f}  "
                        f"{'offset':<10}: {trace_offsets[i]:<12.4f}\n"
                    )
                reference_trace = current_trace.copy()
                reference_trace_num = i + 1
                continue

            corr_ref = np.corrcoef(current_trace, reference_trace)[0, 1]
                       
            if corr_ref < 0:               # negative correlation

                
                offset_jump = ( i >0 and
                        abs(trace_offsets[i] - trace_offsets[i-1]) > 50)

                if offset_jump:

                    reference_trace = current_trace
                    reference_trace_num = i + 1

                    neg_streak = 0        
                    normal.append(all_sac_files[i])
                    normal_offsets.append(trace_offsets[i])
                    normal_corrs.append(abs(corr_ref))

                    f1.write(
                        f"{'trace':<6}: {i+1:<5}  "
                        f"{'reference trace':<16}: {reference_trace_num:<5}  "
                        f"{'corr_normal':<12}: {abs(corr_ref):<8.4f}  "
                        f"{'offset':<10}: {trace_offsets[i]:<12.4f}  "
                        
                    )        

                
                else:
                    
                    neg_streak += 1
                    if neg_streak >= neg_streak_threshold:        # continuous negative → NEW VALID SEGMENT 

                        prev_reference_num = reference_trace_num
                        reference_trace = current_trace.copy()        # treat as new segment instead of removing
                        reference_trace_num = i + 1
                        neg_streak = 0

                        normal.append(all_sac_files[i])
                        normal_offsets.append(trace_offsets[i])
                        normal_corrs.append(corr_ref*(-1))

                        f1.write(
                            f"{'trace':<6}: {i+1:<5}  "
                            f"{'reference trace':<16}: {prev_reference_num:<5}  "
                            f"{'corr_normal':<12}: {corr_ref*(-1):<8.4f}  "
                            f"{'offset':<10}: {trace_offsets[i]:<12.4f}\n"
                        )
                
                    else:
                        
                        polarity_reversed.append(all_sac_files[i])
                        reversed_offsets.append(trace_offsets[i])
                        reversed_corrs.append(corr_ref)

                        f1.write(
                            f"{'trace':<6}: {i+1:<5}  "
                            f"{'reference trace':<16}: {reference_trace_num:<5}  "
                            f"{'corr_normal':<12}: {corr_ref:<8.4f}  "
                            f"{'offset':<10}: {trace_offsets[i]:<12.4f}  "
                            f"REMOVED\n"
                        )

                continue

            
            
            
            neg_streak = 0           #  positive correlation 

            if in_trend or corr_ref >= high_corr_threshold:

                normal.append(all_sac_files[i])
                normal_offsets.append(trace_offsets[i])
                normal_corrs.append(corr_ref)

                f1.write(
                    f"{'trace':<6}: {i+1:<5}  "
                    f"{'reference trace':<16}: {reference_trace_num:<5}  "
                    f"{'corr_normal':<12}: {corr_ref:<8.4f}  "
                    f"{'offset':<10}: {trace_offsets[i]:<12.4f}\n"
                )
                reference_trace = current_trace
                reference_trace_num = i + 1

                
            else:

                polarity_reversed.append(all_sac_files[i])
                reversed_offsets.append(trace_offsets[i])
                reversed_corrs.append(corr_ref)

                f1.write(
                    f"{'trace':<6}: {i+1:<5}  "
                    f"{'reference trace':<16}: {reference_trace_num:<5}  "
                    f"{'corr_normal':<12}: {corr_ref:<8.4f}  "
                    f"{'offset':<10}: {trace_offsets[i]:<12.4f}  "
                    f"REMOVED\n"
                )

    polarity_reversed_array = np.array([sac_file in polarity_reversed for sac_file in all_sac_files])    
    outlier_reverse = polarity_reversed_array
    outlier_reverse_files = [
    sac_file for sac_file, is_outlier in zip(all_sac_files, outlier_reverse)
    if is_outlier]
    
    for sac_file in outlier_reverse_files:
        outlier_file_path = os.path.join(folder, sac_file)
        
        if os.path.exists(outlier_file_path):
            print(f"Deleting reverse traces: {outlier_file_path}")
            os.remove(outlier_file_path)

    sac_files = glob.glob(f"{folder}/*.SAC")

    # convert to filenames only
    sac_filenames = [os.path.basename(f) for f in sac_files]

    remaining_files = [
    f for f in sac_filenames
    if f not in outlier_sac_files and f not in polarity_reversed]
    
    sac_files = glob.glob(f"{folder}/*.SAC")
    fig0, ax0= plt.subplots(figsize=(4.0, 2.3))
    fig1 = plt.figure(figsize=(6.5, 4))
    all_st_3 = Stream()
    for i, sac_file in enumerate(sac_files):

        file_path = os.path.join(folder, sac_file)
        st = read(file_path)
        trace = st[0]
        time = np.arange(0, trace.stats.npts / trace.stats.sampling_rate, trace.stats.delta)
        trace.stats.distance = trace.stats.sac.gcarc * 110  
        # trace.stats.distance = trace.stats.sac.gcarc * 110 *1000 
        # trace_offset = trace.stats.distance /1e3
        all_st_3 += trace 
        
    all_st_3.plot(type='section', fig=fig1, time_down=True, scale=2.5, outfile=save_path4)
    ax1=fig1.axes[0]
    ax1.set_ylim(1000, 200)
    ax1.set_xlabel('Epicentral Distance (Km)', fontsize=10)
    ax1.set_ylabel('Time (s)', fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=9)
    plt.tight_layout()
    plt.savefig(save_path4)
    plt.close(fig1)

    ax0.scatter(normal_offsets, normal_corrs, marker= 'o',s=6, color='black',label='Normal') 
    ax0.scatter(reversed_offsets, reversed_corrs, color='red',s=6, label='Reversed')   
    ax0.set_xlabel("Epicentral Distance (Km)", fontsize=6)
    ax0.set_ylabel("Correlation values", fontsize=6)
    # ax0.set_title("Correlation vs Epicentral Distance", fontsize=10)
    # ax0.legend(fontsize=8, loc="lower left", frameon= False)
    # plt.text(0.02, 1.065, '(a)',transform=plt.gca().transAxes, fontsize=10, va='top', ha='left')
    plt.tick_params(axis='both', which='major', labelsize=5)
    plt.tight_layout()
    plt.savefig(save_path8,format= "pdf",  dpi=300, bbox_inches='tight' )
    plt.close(fig0)
        
    
for current_path in root_path.iterdir():
    if current_path.is_dir():
        for sub_path in current_path.iterdir():
            if sub_path.name == 'LHZ' and sub_path.is_dir():
                target_dir = sub_path/'f-0.0225'
                if target_dir.is_dir():
                    new_subfolder = new_root_path / target_dir.relative_to(root_path)
                    new_subfolder.mkdir(parents=True, exist_ok=True)
                    for sac_file in target_dir.glob('*.sac'):
                        shutil.copy(sac_file, new_subfolder)

                    plot_filename1 = new_subfolder / f'trace_with_envelope.svg'
                    plot_filename2 = new_subfolder / f'observed_&_predicted.svg'
                    scatter_filename = new_subfolder / f'Residual_plot.svg'
                    output_txt_file = new_subfolder / 'Lags.txt'
                    plot_filename3 = new_subfolder/ f'final_traces.svg'
                    #plot_filename4= new_subfolder/ f"removed_traces.png"
                    plot_filename5= new_subfolder/ f"peak_time.pdf"
                    plot_filename6= new_subfolder/ f"trace_with_outlier.svg"
                    plot_filename7= new_subfolder/ f"Correlation Vs Offset.pdf"
                    plot_filename8= new_subfolder/ f"Residual_Trace.svg"
                    plot_filename9= new_subfolder/ f"normal_Trace.svg"
                    save_plots_with_optimized_inversion(new_subfolder, plot_filename1, plot_filename2, scatter_filename, plot_filename3, 
                                                        plot_filename5, plot_filename6, plot_filename7, plot_filename8, plot_filename9)
                    

            elif sub_path.name == 'LHT' and sub_path.is_dir():
                target_dir = sub_path/'f-0.0225'
                if target_dir.is_dir():
                    new_subfolder = new_root_path / target_dir.relative_to(root_path)
                    new_subfolder.mkdir(parents=True, exist_ok=True)

                    for sac_file in target_dir.glob('*.sac'):
                        shutil.copy(sac_file, new_subfolder)
                
                    plot_filename1 = new_subfolder / f'trace_with_envelope.svg'
                    plot_filename2 = new_subfolder / f'observed_&_predicted.svg'
                    scatter_filename = new_subfolder / f'Residual_plot.svg'
                    output_txt_file = new_subfolder / 'lags.txt'
                    plot_filename3 = new_subfolder/ f'final_traces.svg'
                    #plot_filename4= new_inner_subfolder/ f"removed_traces.png"
                    plot_filename5= new_subfolder/ f"peak_time.pdf"
                    plot_filename6= new_subfolder/ f"trace_with_outlier.svg"
                    plot_filename7= new_subfolder/ f"Correlation Vs Offset.pdf"
                    plot_filename8= new_subfolder/ f"Residual_Trace.svg"
                    plot_filename9= new_subfolder/ f"normal_Trace.svg"
                    save_plots_with_optimized_inversion(new_subfolder, plot_filename1, plot_filename2, scatter_filename, plot_filename3, 
                                                                plot_filename5, plot_filename6,plot_filename7,plot_filename8, plot_filename9 )
                    
