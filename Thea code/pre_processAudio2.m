function [timeDomainHRV1, timeDomainHRV2,timeDomainHRV3, timeDomainHRV4] =pre_processAudio2(testnumber, numFiles, folder)

%testnumber = int2str(testnumber);
missingHS = zeros(1,2);

[files yn] = filesForCalculation(testnumber, numFiles);
heartsounds = cell(size(files,1),4);

if yn == 'n'
    return
end
%% figures
areaSegFig = figure('Name','HS After Area Segmentation');
S1_seg = {};
S2_seg = {};
for file = 1:size(files,1)
    %% load data 
    filename = string(files(file,1));
    [~, name, ~] = fileparts(filename)
    disp('Processing: '); disp(name);
    
    audio_resampled = cell2mat(files(file,2));
    resample_fs = cell2mat(files(file,3));
    heartRate = cell2mat(files(file,4));
    sysTimeIntv = cell2mat(files(file,5));
    
    dt2 = 1/resample_fs;
    t2 = 0:dt2:(length(audio_resampled)*dt2)-dt2;
    
    %% removing friction spikes from signal using moving average filter 
    windowLength = dt2*5;
    windowsize = ceil(windowLength/dt2);
    de_spiked_audio = medfilt2(audio_resampled, [windowsize 1]);

    % normalise signal
    de_spiked_audio = de_spiked_audio./max(abs(de_spiked_audio));
    
    despike_fig = figure(20);
    subplot(numFiles,1,file);
    plot(t2,de_spiked_audio); 
    title(string(name))
    hold on;
   
    %% apply envelope thresholding to identify peaks
    [thresh_hilbEnergy,t3] = hilbertEnergyEnvelope2(de_spiked_audio,resample_fs, heartRate);
   
    %% applying area segmentation
    segmentedPeaks = areaSegmentation(thresh_hilbEnergy, resample_fs, t3);
    
    %% identify s1 and s2
    [s1, s2, s1_peaks, s2_peaks, tmp_missingHS] = finds1s2_4(segmentedPeaks,resample_fs,t3,heartRate,sysTimeIntv);
    
    S1_seg{file,1} = name
    S2_seg{file,1} = name
    S1_seg{file,2} = s1(:, 1)
    S2_seg{file,2} = s2(:, 1)
    S1_seg{file,3} = s1(:, 2)
    S2_seg{file,3} = s2(:, 2)
    save(fullfile(folder,'S1_seg.mat'), 'S1_seg')
    save(fullfile(folder,'S2_seg.mat'), 'S2_seg')
    % update number of missing heart sounds counter
    missingHS = missingHS + tmp_missingHS;
  
    % saving s1 and s2 into cells
    heartsounds{file,1} = s1; heartsounds{file,2} = s2;
    heartsounds{file,3} = s1_peaks; heartsounds{file,4} = s2_peaks;
    display(missingHS)

    % Displaying results
    %figure(areaSegFig);
    %subplot(size(files,1),1,file);
    fig = figure();
    plot(t2,audio_resampled./max(abs(audio_resampled))); hold on;
    %plot(t3,thresh_hilbEnergy,'g','LineWidth',1); hold on;
    %plot(t3,segmentedPeaks,'m','LineWidth',1); 
    for i = 1:size(s1,1)
        hold on; xline(s1(i,1),'r','LineWidth',1);
        hold on; xline(s1_peaks(i),'-m');
        hold on; xline(s1(i,2),'r','LineWidth',1);
    end

    for i = 1:size(s2,1)
        hold on; xline(s2(i,1),'g','LineWidth',1);
        hold on; xline(s2_peaks(i),'-y');
        hold on; xline(s2(i,2),'g','LineWidth',1);
    end
    xlabel('Time(s)');
    ylabel('Normalised Amplitude');
    title(name); % make own legend
    f=gca;
    exportgraphics(f, fullfile(folder, strcat(string(name), "_seg.png")));
end
exportgraphics(despike_fig, fullfile(folder, strcat(name, "-despiked.png")));

% calculate HRV metrics
disp('------------------------'); 
[timeDomainHRV1, timeDomainHRV2] = calculateTimeDomainHRVSounds(heartsounds,missingHS);
disp('------------------------'); 
[timeDomainHRV3, timeDomainHRV4] = calculateTimeDomainHRVPeaks(heartsounds,missingHS);
disp('------------------------'); 

% exporting data into excel
% filename = 'HRV.xlsx';
% writematrix(timeDomainHRV1,filename,'Sheet',1,'Range','E1:J1');
% writematrix(timeDomainHRV2,filename,'Sheet',2,'Range','E1:J1');
% writematrix(timeDomainHRV3,filename,'Sheet',3,'Range','E1:J1');
% writematrix(timeDomainHRV4,filename,'Sheet',4,'Range','E1:J1');

end
