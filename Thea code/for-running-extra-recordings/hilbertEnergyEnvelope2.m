function [thresh_hilbEnergy, t] = hilbertEnergyEnvelope(audio_signal,fs,hr)

% method is based on paper by Sharma et al (referring to Sharma et al 'An
% Algorithm for Heart Rate Extraction from Acoustic Recordings at the
% Neck')


%% for testing function
% [audio_signal fs hr] =pre_processAudio2(["relaxed-1.wav"],3)

%% figures

% hilbGraph = figure('Name', 'Hilbert Energy Envelopes');
% hilbGraphPostThresh1 = figure('Name', 'Hilbert Envelope After Initial Threshold');
% hilbGraphPostThresh2 = figure('Name', 'Hilbert Energy Envelope Post Adaptive Threshold');


%% splitting data into blocks of 2s (can we change this to suit the length of the recording?)

% get new time values
dt = 1/fs;
tot_t = (length(audio_signal)*dt)-dt;
t = 0:dt:tot_t;

% take hilbert transform
y = hilbert(audio_signal);

% calculate hilbert envelope
inst_amp = sqrt(audio_signal.^2 + y.^2);
hilbEnergy = abs(inst_amp.^2);

% applying 1D moving average filter
windowSize = 50; 
b = (1/windowSize)*ones(1,windowSize); a = 1;
filt_hilbEnergy = filter(b,a,hilbEnergy);

% pad the end with zeros for 7s in total
pad_zeros = zeros(7000-length(filt_hilbEnergy),1);
filt_hilbEnergy = [filt_hilbEnergy; pad_zeros];
t = [t t(end)+dt:dt:7-dt];

% separate envelope and time array into 2s intervals which overlap by 50%
filt_hilbEnergy = buffer(filt_hilbEnergy,2000,1000);
filt_hilbEnergy(:,1) = [];
t = buffer(t',2000,1000);
t(:,1) = [];

thresh_hilbEnergy = zeros(size(filt_hilbEnergy));

% number of segments:
L = size(t,2);

%% For each segment: ---------------------------------------------------------------------
thresh_arr = zeros(L,1);
thresh_idx = zeros(size(t));

% uper and lower bounds determined from heart rate approximated by autocorrelation
lbound = floor((hr*4)/60)+1; % higher boundaries to ensure smaller beats aren't excluded
hbound = ceil((hr*4)/60)+4;


for z = 1:L
    % initial threshold 
    thresh = 0.1; 

    % taking normalised amplitude values
    filt_hilbEnergy(:,z) = filt_hilbEnergy(:,z)./max(filt_hilbEnergy(:,z));
    thresh_hilbEnergy(:,z) = filt_hilbEnergy(:,z);

    %% plot hilbert transform for each segment ----------------------------------- 
%     figure(hilbGraph);
%     subplot(1,L,z);
%     plot(t(:,z), filt_hilbEnergy(:,z));
%     xlabel('Time');ylabel('Amplitude');
%     title(['Segment:',num2str(z)]);
    
    %% applying intial threshold ---------------------------------------------------
    
    T = thresh*max(filt_hilbEnergy(:,z));
    thresh_idx = find(filt_hilbEnergy(:,z) < T);
    thresh_hilbEnergy(thresh_idx, z) = 0;
    
    %% plot initial 0.20 thresholded --------------------------------------------
    
%     figure(hilbGraphPostThresh1);hold on;
%     subplot(1,L,z);
%     plot(t(:,z), thresh_hilbEnergy(:,z));hold on;
%     xlabel('Time');ylabel('Amplitude');
%     title(['Segment ',num2str(z), ' & Thresh ', num2str(thresh)]);


%% applying adaptive thresholding ----------------------------------------------------
    
    peaks_z = countPeaks(thresh_hilbEnergy(:,z));
    
    iter_count = 1;
    while ((peaks_z < lbound) || (peaks_z > hbound))

        % change t according to number of peaks§
        if ((peaks_z < lbound) && (thresh > 0))
            thresh = thresh - 0.02;
        else
            if thresh <= 0.5
                thresh = thresh + 0.02;
            end
        end
        
        if (round(thresh,3) <= 0) || (thresh >= 0.5)
            thresh = 0.01;
            break;
        end

        % re-starting thresh_hilbEnergy for that segment
        thresh_hilbEnergy(:,z) = filt_hilbEnergy(:,z); 
        
        % apply threshold
        T = thresh*max(filt_hilbEnergy(:,z));
        thresh_idx = find(filt_hilbEnergy(:,z) < T);
        thresh_hilbEnergy(thresh_idx, z) = 0;
        
        % store new threshold in array
        thresh_arr(z) = thresh;
        
        % re-calculate number of peaks
        peaks_z = countPeaks(thresh_hilbEnergy(:,z));
        
        % update boundaries if number of iterations becomes too great
        iter_count = iter_count+1;
        if iter_count > 10
            lbound = lbound - 1;
            break;
        end
    end
    
    % storing final threshold into array
    thresh_arr(z) = thresh;

    
    %% plot thresholded Hilbert Energy Envelope
%     figure(hilbGraphPostThresh2);
%     subplot(1,L,z);
%     plot(t(:,z), thresh_hilbEnergy(:,z));hold on;
%     xlabel('Time');ylabel('Amplitude');
%     title(['Segment ',num2str(z), ' & Thresh ', num2str(thresh_arr(z))]);
%     
    
end


%% displaying results ---------------------------------------------------------------------------

% plot original signal against hilbert envelope
% figure(hilbGraphPostThresh2); subplot(2,1,1);
% plot([0:dt:(length(audio_signal)*dt)-dt],audio_signal./max(abs(audio_signal)),'y'); hold on;
% plot(t,filt_hilbEnergy,'k', 'LineWidth', 1);
% title('Hilbert Envelope');
% 
% % plot original signal against thresholded envelope
% hold on; subplot(2,1,2);
% plot([0:dt:(length(audio_signal)*dt)-dt],audio_signal./max(abs(audio_signal))); hold on;
% plot(t,thresh_hilbEnergy,'k', 'LineWidth', 1);
% xlabel('Time(s)');
% ylabel('Amplitude');
% title('Hilbert Transform ');



end


%% ---------------------------------------------------------------------------------------


function [numPeaks] = countPeaks(signal)
%% counts the number of peaks in a thresholded signal
    ne0 = find(signal~=0);                                   % Nonzero Elements
    ix0 = unique([ne0(1); ne0(diff([0; ne0])>1)]);  

    numPeaks = length(ix0);

end