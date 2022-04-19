%% Automated function to test .wav files for Aortic Stenosis
function [final] = finalFunction(myFolder,op)

% Get a list of all files in the folder with the desired file type
filePattern = fullfile(myFolder, '*.wav');
Files = dir(filePattern);
for k = 1:length(Files)
    baseFileName = Files(k).name
    fullFileName = fullfile(Files(k).folder, baseFileName);

    %% feeds file to be filtered and removes murmur
    [orig,filt,f_data,fs,total_t] = murmurRemoval(fullFileName,150,250);
    % Returns:
    %   1). orig = original recording
    %   2). filt = filtered data through 150-250Hz bandpass
    %   3). f_data = filtered and murmur removed data for boundary detection
    
    %% Iterate to find the best threshold for the findPeaks function
    ratioNew = [0.1, 0.15, 0.2, 0.25, 0.3];
    for l = 1:length(ratioNew)
        rat = ratioNew(l);
        
        % Return initial boundaries and peaks
        [heart_sounds,check] = hilbertAndFindPeaks(f_data,fs,50,0.1,rat,0,0,total_t);
        % heart_sounds = matrix of boundary start time, heart sound classification, boundary end time for each heart sound

        % if hilbertAndFindPeaks can't find boundaries, script jumps to next file
        if check == 0
            continue
        end

        %% Accurate hearbeats based off of time intervals between S1 & S2
        % Calculate the number of accurate heartbeats
        [mean_t,number(k,1),s1_m(k,1),s1_std(k,1),s2_m(k,1),s2_std(k,1)] = AcceptableWindows(heart_sounds);
        % Returns variables:
        %   1). mean_t = mean heartbeat length in seconds of acceptable heartbeats
        %   2). number = number of acceptable heartbeats
        %   3). s1_std & s2_std = number of acceptable heartbeats
        
        %% Calculates amplitudes between heart sounds & provides locations of boundaries
        [ampt_12,ampt_21,S12_t,S21_t,ratio] = CalculateAmplitude(heart_sounds,fs,orig,filt);
        % Returns variables:
        %   1). ampt_12 & ampt_21 = systolic and diastolic amplitudes
        %   2). S12_t & S21_t = column matrix of the times at which the peaks occur
        %   3). ratio = colum matrix of the systolic to diastolic ratio of each individual cycle

        %% Calculate heartbeat metrics
        % number of diastole/systole visible
        num_s = length(ampt_12); num_d = length(ampt_21);  
        % number of heartbeats visible 
        if num_s > num_d
            num(k,1) = num_d;
        else
            num(k,1) = num_s;
        end

        % percentage of acceptable heartbeats from visible ones
        pc(k,l) = (number(k,1)/num(k,1))*100;
    end
    
    %% Calculate best threshold ratio to use for findPeaks function
    bestR = zeros(1);
    
    % find index of best threshold
    val = max(pc(k,:)); idx = find(pc(k,:)==val);
    pcn(k) = val;
    
    % bestR is the optimised threshold
    bestR = ratioNew(idx(1));
    
    %% Run whole algorithm again using best threshold for optimum results
    [heart_sounds,check] = hilbertAndFindPeaks(f_data,fs,50,0.1,bestR,0,0,total_t);
    
    % if hilbertAndFindPeaks can't find boundaries, script jumps to next file
    if check == 0
            continue
    end
    [ampt_12,ampt_21,S12_t,S21_t,ratio] = CalculateAmplitude(heart_sounds,fs,orig,filt);
    
    %% Expected number of heart beats based off of length of recording
    % hcl = expected heart cycle length = 0.9 seconds
    hcl = 0.9;
    num_rough(k,1) = floor(total_t/hcl);
    
    %% Calculate amplitude metrics
    [sd_ratio(k,1),avg_sys_amp(k,1),avg_dias_amp(k,1),~] = AmplitudeMetrics(ampt_12,ampt_21,ratio);
    % Returns:
    %   1). sd_ratio = systolic to diastolic ratio
    %   2). avg_sys_amp = average systolic amplitude
    %   3). avg_dias_amp = average diastolic amplitude

    %% filename
    B(k,:) = baseFileName(1:5);
    B(k,:) = convertCharsToStrings(B(k,:));
    
end

%% Print variables to a table
varNames = ["File", "S/D A Ratio", "Expected HB","HB Seen", "Acceptable HB", "S1 int Mean","S2 int Mean"];
T = table(B,round(sd_ratio,2),num_rough,num,number,round(s1_m,2),round(s2_m,2), 'VariableNames', varNames);

%% final results & metrics
new_sd_ratio = round(sd_ratio,2);

% finding all sensitivities and specifities for each threshold
if op == 1
    x = 1;
    for p = 0.8:0.05:2.85
        total = new_sd_ratio(new_sd_ratio<p);
        final(x) = length(total)/length(new_sd_ratio);
        x = x+1;
    end 
elseif op == 2
    y = 1;
    for p = 0.8:0.05:2.85
        total = new_sd_ratio(new_sd_ratio>p);
        final(y) = length(total)/length(new_sd_ratio);
        y = y+1;
    end 
end
