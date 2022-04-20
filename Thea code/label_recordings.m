clc
clearvars
disp_ECG = 1
label_data = 0

PCG_dir = "../../EKO/Patient 1 data/PCG";
PCG_files = dir(fullfile(PCG_dir, '*.wav'));
[~,p_ind]=sort({PCG_files.name});
PCG_files = PCG_files(p_ind);

if disp_ECG 
    ECG_dir = "../../EKO/Patient 1 data/ECG";
    ECG_files = dir(fullfile(ECG_dir, '*.wav'));
    [~,e_ind]=sort({ECG_files.name});
    ECG_files = ECG_files(e_ind);
end


for i = 1:size(PCG_files,1)
  PCG_paths(i) = fullfile(PCG_dir, PCG_files(i).name); % store filenames in str
  if disp_ECG 
    ECG_paths(i) = fullfile(ECG_dir, ECG_files(i).name);
  end
end

%% Global Variables
S1 = {};
S2 = {};

%% run through recordings for manual annotation
for i=1:length(PCG_paths)
    figure();
    disp(PCG_files(i).name);
    resample_rate = 1000;
    % read audio file 
    [p_signal, p_fs] = audioread(PCG_paths(i));

    % re-sample file
    p_signal = resample(p_signal, resample_rate, p_fs);
    dt = 1/resample_rate; 
    t = 0:dt:(length(p_signal)*dt)-dt;

    % de-noise signal
    p_signal = applyButterworthBandpassFilter(25, 165, 3, resample_rate, p_signal);
    p_signal = p_signal./max(abs(p_signal));
    
    % display hilbert envelope to see if there are distinct peaks in env
    hilb = hilbert_envelope(p_signal);
     
    if disp_ECG
        [e_signal, e_fs] = audioread(ECG_paths(i));
        e_dt = 1/e_fs; 
        e_t = 0:e_dt:(length(e_signal)*e_dt)-e_dt;

        e_signal = applyButterworthBandpassFilter(1.5, 100, 4, e_fs, e_signal);
        
        wo = 50/(e_fs/2);  
        bw = wo;
        [num,den] = iirnotch(wo,bw);
        e_signal = filtfilt(num, den, e_signal);
        
        windowWidth = 10; % Whatever you want.
        kernel = ones(windowWidth,1) / windowWidth;
        e_signal = filtfilt(kernel, 1, e_signal);

        e_signal = e_signal./max(abs(e_signal));

        % display cleaned signal 
        fig = tiledlayout(3,1);
        fig.TileSpacing = 'compact';
        fig.Padding = 'compact';
        title(fig, PCG_files(i).name)
        ax1 = nexttile([1 1]);
        plot(e_t, e_signal);
        ax2 = nexttile([2 1]);
        plot(t, p_signal, t, hilb);
    else
        fig = figure();
        title(PCG_files(i).name);
        plot(t, p_signal);
    end

    if label_data
        good_rec = input(['Is recording good enough?'], 's')
        if good_rec == 'y'
            S1{i,1} = PCG_files(i).name;
            S2{i,1} = PCG_files(i).name;
    
            num = input('How many S1 sounds do you see? ');
            disp('Choose start of S1 sounds');
            [s1_start,y] = ginput(num);
            % store value of x coordinate
            S1{i,2} = s1_start;
            
            %%
            disp('Choose end of S1 sounds');
            [s1_end,y] = ginput(length(s1_start));
            % store value of x coordinate
            S1{i,3} = s1_end;
            
            %%
            num = input('How many S2 sounds do you see? ');
            disp('Choose start of s2 sounds');
            [s2_start,y] = ginput(num);
            % store value of x coordinate
            S2{i,2} = s2_start;
            
            %%
            disp('Choose end of s2 sounds');
            [s2_end,y] = ginput(length(s2_start));
            % store value of x coordinate
            S2{i,3} = s2_end;
        
            save(fullfile(PCG_dir,'S1.mat'), 'S1')
            save(fullfile(PCG_dir,'S2.mat'), 'S2')
            
        
            y = input(['Do you wish to continue to the next recording? '],'s');
            if y == 'y'
                continue;
            else
                disp(['Finished at file: ', num2str(i)]);
                break;
            end
        end
    end
end

function hilbEnergy = hilbert_envelope(signal)

y = hilbert(signal);
inst_amp = sqrt(signal.^2 + y.^2);
hilbEnergy = abs(inst_amp.^2);
windowSize = 50; 
b = (1/windowSize)*ones(1,windowSize); a = 1;
hilbEnergy = filtfilt(b,a,hilbEnergy);

end