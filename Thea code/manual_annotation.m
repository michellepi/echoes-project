clear;clc;
% script for manually annotating recordings

%% load filenames
Folder = "../../sound_files/good-recordings"; % add path to folder containing .wav files here


% get filenames for signal extraction 
FileList = dir(fullfile(Folder, '*.wav'));
filename = strings(size(FileList,1),1); % string array to store the filenames of the recordings to be segmented
for i = 1:size(FileList,1)
  filename(i) = fullfile(Folder, FileList(i).name); % store filenames in str
end


%% Global Variables
fs2 = 1000;
S1 = cell(length(filename),2);
S2 = cell(length(filename),2);

%% run through recordings for manual annotation
for i=2:length(filename)
    disp(filename(i));
    % read audio file 
    [signal, fs] = audioread(filename(i));
    % re-sample file
    signal = resample(signal,fs2,fs);
    dt = 1/fs2; t = 0:dt:(length(signal)*dt)-dt;

    % de-noise signal
    signal = applyButterworthBandpassFilter(25, 165, 3, fs2, signal);
    signal = signal./max(abs(signal));
    
    % display hilbert envelope to see if there are distinct peaks in env
    hilb = hilbert_envelope(signal);

    % display cleaned signal 
    fig1 = figure('Name',filename(i));
    plot(t,signal,t,hilb);
    
    %%
    num = input('How many S1 sounds do you see? ');
    disp('Choose start of S1 sounds');
    [s1_start,y] = ginput(num);
    % store value of x coordinate
    S1{i,1} = s1_start;
    
    %%
    disp('Choose end of S1 sounds');
    [s1_end,y] = ginput(length(s1_start));
    % store value of x coordinate
    S1{i,2} = s1_end;
    
    %%
    num = input('How many S2 sounds do you see? ');
    disp('Choose start of s2 sounds');
    [s2_start,y] = ginput(num);
    % store value of x coordinate
    S2{i,1} = s2_start;
    
    %%
    disp('Choose end of s2 sounds');
    [s2_end,y] = ginput(length(s2_start));
    % store value of x coordinate
    S2{i,2} = s2_end;
    
    %%
    y = input('Do you wish to continue to the next recording? ','s');
    if y == 'y'
        continue;
    else
        disp(['Finished at file: ', num2str(i)]);
        break;
    end
end

%% for viewing labels
for n=1:length(S1)
s1 = cell2mat(S1(n,:));
    for i = 1:size(s1,1)
        hold on; xline(s1(i,1),'r','LineWidth',1);
        hold on; xline(s1(i,2),'r','LineWidth',1);
    end
end

for n=1:length(S2)
s2 = cell2mat(S2(n,:));
    for i = 1:size(s2,1)
        hold on; xline(s2(i,1),'g','LineWidth',1);
        hold on; xline(s2(i,2),'g','LineWidth',1);
    end
end





function hilbEnergy = hilbert_envelope(signal)

y = hilbert(signal);
inst_amp = sqrt(signal.^2 + y.^2);
hilbEnergy = abs(inst_amp.^2);
windowSize = 50; 
b = (1/windowSize)*ones(1,windowSize); a = 1;
hilbEnergy = filter(b,a,hilbEnergy);

end