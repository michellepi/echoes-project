function filtered_signal = applyButterworthBandpassFilter(lowbound, highbound, order, sampling_freq, signal)
% Function returns a filtered signal after applying a n-th order bandpass
% butterworth filter for frequencies between lowbound to highbound in the
% unfiltered signal.

% Apply lowpass butterworth filter
[b_low,a_low] = butter(order,2*highbound/sampling_freq,'low');
lowpass_filt = filtfilt(b_low,a_low,signal);

% Apply highpass butterworth filter
[b_high,a_high] = butter(order,2*lowbound/sampling_freq,'high');
filtered_signal = filtfilt(b_high,a_high,lowpass_filt);
end 