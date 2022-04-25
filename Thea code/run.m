%% instructions for running code

% pre_processAudio2 takes one input: filename which is the directory of a
% given .wav file
%   - pre_processAudio2 will run te entire algorithm then produce a figure
%   of the recording with the identified peaks; if you would like to see
%   the peaks after peak selection you can uncomment lines 60 and 61 in
%   pre_processAudio2


% WARNING: the overall code still has major bugs that need fixing and is unable to
% deal with lots of artefacts so it may not work in some cases


% assuming all .wav files are in the same folder, just run a for loop for
% each file (the code below should do the job, run the script to execute for all files in the folder or get rid of
% for loop to run for specific file)

%% ENTER DIRECTORY OF FOLDER CONTAINING RECORDINGS
Folder = ['../../EKO/Patient 2/PCG']; % change this to the full directory of the folder containing all the .wav files

%% creating a matrix to store filenames
FileList = dir(fullfile(Folder, '**', '*.wav'));
display(FileList)
filenames = strings(size(FileList,1),1);
for i = 1:size(FileList,1)
    
  filenames(i) = string(fullfile(FileList(i).folder, FileList(i).name));
end

%% run heart sound segmentation for all files in filenames
%for i=1:size(filenames,1)
 %   pre_processAudio2(filenames(i), 1); 
%end
pre_processAudio2(filenames, size(filenames,1), Folder); % Change second input parameter to set reqquired number of files