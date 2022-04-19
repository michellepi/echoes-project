function [freqDomainHRV1, freqDomainHRV2] = calculateFreqDomainHRV(heartsounds,fs)

numFiles = size(heartsounds,1);

%% calculate HRV paramters for each file
for file = 1:numFiles
    % load data
    s1 = cell2mat(heartsounds(file,1));
    s2 = cell2mat(heartsounds(file,2));
    
    % calculate intervals
    ibi_s1 = diff(s1(:,1)); ibi_s2 = diff(s2(:,1));
    
    % using periodogram to get psd estimate
    psdest = psd(spectrum.periodogram,ibi_s1,'Fs',fs,'NFFT',length(ibi_s1));
    LFpower = avgpower(psdest,[0.04 0.15]);
    HFpower = avgpower(psdest,[0.4 1]);

    xdft = fft(ibi_s1);
    xdft = xdft(1:length(ibi_s1)/2+1);
    xdft(2:end-1) = 2*xdft(2:end-1);
    psdest = 1/(length(ibi_s1)*fs)*abs(xdft).^2;
    freq = 0:fs/length(ibi_s1):fs/2;
    plot(freq/100,10*log10(psdest));
    LF1 = psdest(find(freq> 0.04 & freq<0.15));
    
    figure(14)
    subplot(numFiles,1,numFiles);
    plot(freq/100,10*log10(psdest));
    xlabel('Frequency (Hz)'); ylabel('Power Spectral Density (dB/Hz)');
    
end


end