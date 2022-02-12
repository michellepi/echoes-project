function [required_files yn] = filesForCalculation(testnumber, numFiles)
% function identifies the appropriate files that fall below a noise
% threshold. These files will then be used for the calculation of HRV
% outputs:
%   - required_files: cell array which stores filename,audio_resampled,t2,heart rate and systolice time interval
%   - yn: yes or no from user to indicate continuation in the case of less heart ounds found than orginally asked

%% Read all recordings in folder to assess their noise quality
testnumber = num2str(testnumber); % folder containing recordings
Folder = append('../../sound_files/hs-with-ground-truth-hrv/',testnumber);
FileList = dir(fullfile(Folder, '*.wav'));
filenames = strings(size(FileList,1),1);
for i = 1:size(FileList,1)
  filenames(i) = fullfile(Folder, FileList(i).name);
end

%% global variables
required_files = cell(size(FileList,1),5);
cell_row = 1;

%% for testing 
% append the following in function definition for use in other functions
% for testing: [segmentedPeaks resample_fs t3 audio_resampled t2 sysTimeIntv heartRate] =
% filenames = append('../../sound_files/hs-with-ground-truth-hrv/',testnumber,'13-04-2021_09.43.04.wav');
% filenames = '../../sound_files/hs-with-ground-truth-hrv/1/08-03-2021_16.37.03.wav'
% filenames = '../../sound_files/hs-with-ground-truth-hrv/1/08-03-2021_16.38.58.wav'
% filenames = string(filenames);


for file = 1:length(filenames)
    %% Reading audio files
    disp(filenames(file))
    [audio_file, Fs] = audioread(filenames(file));
    audio_file = audio_file(:,1);

    %% Pre-processing -------------------------------------------------------------
    % resampling to 1000Hz
    resample_fs = 1000; % downsampling reduces the number of samples to be processed at later stages
    audio_resampled = resample(audio_file, resample_fs, Fs);
    
    %% get estimate of heart rate and systolic time interval
    [heartRate, sysTimeIntv] = getHeartRateSchmidt(audio_resampled, resample_fs, false); % what does the given graph mean?
    
    %% apply bandpass filter
    % filter signal to 25 and 165Hz (general range of heart sounds)
    audio_resampled = applyButterworthBandpassFilter(30, 125, 3, resample_fs, audio_resampled);
    
    % apply adaptive thresholding according to Jain et al
    [noise_level, noise_score] = calculateNoiseLevel(audio_resampled);
    disp(['------------------Noise level: ',noise_level]);
    
    if noise_score < 120 && noise_score > 60
        % if recording meets threshold for noise, store it in output array
        required_files{cell_row,1} = filenames(file);
        required_files{cell_row,2} = audio_resampled;
        required_files{cell_row,3} = resample_fs;
        required_files{cell_row,4} = heartRate;
        required_files{cell_row,5} = sysTimeIntv;
        
        cell_row = cell_row + 1;
    else
        continue;
    end


end
if cell_row <= size(required_files,1)
    required_files(cell_row:end,:) = [];
end

if size(required_files,1) > numFiles
    % we must reduce the number of recordings by num_rm
    num_rm = size(required_files,1) - numFiles;
    % choose random permutation of numbers between 1 and size(required_files,1); (always keep first and last file)
    Index = 1+randperm(size(required_files,1)-2, num_rm);
    required_files(Index,:) = [];
    yn = 'y';
elseif size(required_files,1) < numFiles
    % return error message and ask if user still wants to continue
    yn = input(['There are only ',num2str(size(required_files,1)),' files found. Would you like to continue? (y/n) '],'s');
elseif size(required_files,1) == numFiles
    yn = 'y';
end
end

