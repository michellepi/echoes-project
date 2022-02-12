function [timeDomainHRV1, timeDomainHRV2] = calculateTimeDomainHRVPeaks(heartsounds,missingHS)
% function calculates the time-domain parameters of HRV 
% input: 
%   - heartsounds: cell array containing start and end indices of s1 and s2
%   sounds as well as the locations of s1 and s2 peaks from a given number of files
% output: 
%   - timeDomainHRV: time-domain HRV parameters from files

%% for testing:
% [heartsounds, missingHS] = pre_processAudio2(3,5);

%% initialising variables
sum_rmssd1 = 0; sum_rmssd2 = 0;
sum_sdnn1 = 0; sum_sdnn2 = 0;
sum_means1s1 = 0; sum_means2s2 = 0;
sum_pnn501 = 0; sum_pnn502 = 0;
sum_hr1 = 0; sum_hr2 = 0;
size_org1 = 0; size_org2 = 0;
size_corrected1 = 0; size_corrected2 = 0;

numFiles = size(heartsounds,1);

%% calculate HRV paramters for each file
for file = 1:numFiles
    % load data
    s1 = cell2mat(heartsounds(file,3)); size_org1 = size_org1 + size(s1,1);
    s2 = cell2mat(heartsounds(file,4)); size_org2 = size_org2 + size(s2,1);
    
    % calculate inter-beat intervals (IBI) for each heart sound
    ibi_s1 = diff(s1(:)); ibi_s2 = diff(s2(:));
    
    % removing outliers
    ibi_s1 = rmoutliers(ibi_s1); ibi_s2 = rmoutliers(ibi_s2);
    ibi_s1(ibi_s1>1.2) = []; ibi_s1(ibi_s1<0.5) = [];
    ibi_s2(ibi_s2>1.2) = []; ibi_s2(ibi_s2<0.5) = [];
%     disp(ibi_s1);disp(ibi_s2);
    size_corrected1 = size_corrected1 + length(ibi_s1); size_corrected2 = size_corrected2 + length(ibi_s2);
    
    % calculate differences in IBI 
    diff_ibi1 = abs(diff(ibi_s1)); diff_ibi2 = abs(diff(ibi_s2));
    
    % calculate rmssd of file and sum
    rmssd_s1 = rms(diff_ibi1); rmssd_s2 = rms(diff_ibi2);
    sum_rmssd1 = sum_rmssd1 + rmssd_s1;
    sum_rmssd2 = sum_rmssd2 + rmssd_s2;
    
    % calculate sdnn of file and sum
    sdnn_s1 = std(ibi_s1); sdnn_s2 = std(ibi_s2);
    sum_sdnn1 = sum_sdnn1 + sdnn_s1;
    sum_sdnn2 = sum_sdnn2 + sdnn_s2;
    
    % find mean s1s1 and s2s2 interval of file and take sum
    mean_s1s1 = mean(ibi_s1); mean_s2s2 = mean(ibi_s2);
    sum_means1s1 = sum_means1s1 + mean_s1s1;
    sum_means2s2 = sum_means2s2 + mean_s2s2;
    
    % calculate pNN50 and sum
    pnn501 = size(find(ibi_s1 > 0.050),1)/size(s1,1);
    pnn502 = size(find(ibi_s2 > 0.050),1)/size(s2,1);
    sum_pnn501 = sum_pnn501 + pnn501;
    sum_pnn502 = sum_pnn502 + pnn502;
    
    % find hr from file
    hr_s1 = 60/mean_s1s1; hr_s2 = 60/mean_s2s2;
    sum_hr1 = sum_hr1 + hr_s1;
    sum_hr2 = sum_hr2 + hr_s2;
      
end 

%% calculate overall time-domain parameters 
% calculate mean RMSSD
RMSSD1 = sum_rmssd1/numFiles; RMSSD2 = sum_rmssd2/numFiles;

% calculate mean SDNN
SDNN1 = sum_sdnn1/numFiles; SDNN2 = sum_sdnn2/numFiles; 

% calculate mean S1S1 and S2S2 intervals
MEANS1S1 = sum_means1s1/numFiles; MEANS2S2 = sum_means2s2/numFiles;

% calculate mean pNN50
PNN501 = sum_pnn501/numFiles; PNN502 = sum_pnn502/numFiles;

% calculate mean heart rate
HR1 = sum_hr1/numFiles; HR2 = sum_hr2/numFiles;

% storing for output
timeDomainHRV1 = [HR1, RMSSD1*1000, SDNN1*1000, PNN501*100, MEANS1S1*1000, missingHS(1)+(size_org1-size_corrected1)];
timeDomainHRV2 = [HR2, RMSSD2*1000, SDNN2*1000, PNN502*100, MEANS2S2*1000, missingHS(1)+(size_org2-size_corrected2)];

%% display values
% disp(['Mean RMSSD (from s1): ', num2str(RMSSD1*1000),' ms']);
% disp(['Mean RMSSD (from s2): ', num2str(RMSSD2*1000),' ms']);
% 
% disp(['Mean SDNN (from s1): ',num2str(SDNN1*1000),' ms']);
% disp(['Mean SDNN (from s2): ',num2str(SDNN2*1000),' ms']);
% 
% disp(['Mean S1S1 interval: ',num2str(MEANS1S1*1000),' ms']);
% disp(['Mean S1S1 interval: ',num2str(MEANS2S2*1000),' ms']);
% 
% disp(['Mean PNN50 (from s1): ',num2str(PNN501*100),'%']);
% disp(['Mean PNN50 (from s2): ',num2str(PNN502*100),'%']);
% 
% disp(['Mean HR (S1)', num2str(HR1),'bpm']);
% disp(['Mean HR (S2)', num2str(HR2),'bpm']);
% 
% disp(['Number of missing S1 sounds: ',num2str(missingHS(1))]);
% disp(['Number of missing S2 sounds: ',num2str(missingHS(2))]);

end