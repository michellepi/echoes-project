clear;clc;
% script for viewing recordings
%% load filenames
Folder_ecg = "../../Patient 1 data/ECG"; % add path to folder containing .wav files here
Folder_pcg = "../../Patient 1 data/PCG"
% get filenames for signal extraction 
FileList = dir(fullfile(Folder, '*.wav'));
filename = strings(size(FileList,1),1); % string array to store the filenames of the recordings to be segmented
for i = 1:size(FileList,1)
  filename(i) = fullfile(Folder, FileList(i).name); % store filenames in str
end

%%
clean_as_filenames = strings();
fs2 = 1000;
idx_file = 1;

% uncomment to view labels for good recordings
%load('good-as-rec-labels.mat'); % Contains S1 and S2 cells with annotations of onset and offset of HSs


%% view each recording
for i = 1:length(filename)
    disp(filename(i));
    % read audio file 
    [signal, fs] = audioread(filename(i));
    % re-sample file
    signal = resample(signal,fs2,fs);
    dt = 1/fs2; t = 0:dt:(length(signal)*dt)-dt;

    % de-noise signal
    signal = applyButterworthBandpassFilter(25, 165, 3, fs2, signal);
    signal = signal./max(abs(signal));
    
    hilb = hilbert_envelope(signal);

    % display cleaned signal 
    fig1 = figure('Name',filename(i));
    plot(t,signal,t,hilb);
    
    % ask user if file is clean enough for as detection
%     y = input('Is this clean enough? ','s');
%     if y == 'y'
%         clean_as_filenames(idx_file) = filename(i);
%     end 
%     idx_file = idx_file + 1;

    
    % uncomment to view labels 
   % s1 = cell2mat(S1(i,:));
    %for n = 1:size(s1,1)
        %hold on; xline(s1(n,1),'r','LineWidth',1);
        %hold on; xline(s1(n,2),'r','LineWidth',1);
   % end

    %s2 = cell2mat(S2(i,:));
    %for n = 1:size(s2,1)
        %hold on; xline(s2(n,1),'g','LineWidth',1);
        %hold on; xline(s2(n,2),'g','LineWidth',1);
    %end
    
end

function hilbEnergy = hilbert_envelope(signal)

y = hilbert(signal);
inst_amp = sqrt(signal.^2 + y.^2);
hilbEnergy = abs(inst_amp.^2);
windowSize = 50; 
b = (1/windowSize)*ones(1,windowSize); a = 1;
hilbEnergy = filter(b,a,hilbEnergy);

end