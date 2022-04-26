clc
clearvars
disp_ECG = 1

label_data = 0
root = "../../Healthy Volunteer/EMILIE - SITTING";

PCG_dir = fullfile(root, "/PCG")
PCG_files = dir(fullfile(PCG_dir, '**', '*.wav'));
[~,p_ind]=sort({PCG_files.name});
PCG_files = PCG_files(p_ind);

if disp_ECG 
    ECG_dir = fullfile(root, "/ECG")
    ECG_files = dir(fullfile(ECG_dir,'**', '*.wav'));
    [~,e_ind] = sort({ECG_files.name});
    ECG_files = ECG_files(e_ind);
end

for i = 1:size(PCG_files,1)
  file = fullfile(PCG_files(i).folder, PCG_files(i).name)
  PCG_paths(i) = string(file); % store filenames in str
  
  if disp_ECG 
      file = fullfile(ECG_files(i).folder, ECG_files(i).name)
    ECG_paths(i) = string(file);
  end
end

%% Global Variables
S1 = {};
S2 = {};

%% run through recordings for manual annotation
for i=1:length(PCG_paths)
    disp(PCG_files(i).name);
    resample_rate = 1000;
    % read audio file 
    disp(PCG_paths(i))
    [p_signal, p_fs] = audioread(PCG_paths(i));

    M=p_fs;% define the window-length
    window=hanning(M); % define the window to use (Hanning)
    overlap=0.5; % define the overlap for windows
    [P,f]=spec2(p_signal,window,overlap,p_fs);% calculate the PSD
    figure();
    plot(f,P);% plot the result
    xticks(0:100:max(f));
    xlabel('frequency (Hz)');% label axes
    ylabel('PSD (microV^2/Hz)');
    p_title = strcat('PSD of PCG signal ', PCG_files(i).name)

    title(p_title);

    % re-sample file
    p_signal = resample(p_signal, resample_rate, p_fs);
    dt = 1/resample_rate; 
    t = 0:dt:(length(p_signal)*dt)-dt;

    % de-noise signal
    %p_signal = applyButterworthBandpassFilter(25, 165, 3, resample_rate, p_signal);
     p_low={};
     p_low.N=2;
     p_low.type='butter';
     p_low.lphp='low';
     [b,a]=make_digital_filter(150, p_fs, p_low);
     p_signal = filtfilt(b, a, p_signal);

     p_high={};
     p_high.N=2;
     p_high.type='butter';
     p_high.lphp='high';
     [b,a]=make_digital_filter(25, p_fs, p_high);
     p_signal = filtfilt(b, a, p_signal);

    p_signal = p_signal./max(abs(p_signal));
    
    % display hilbert envelope to see if there are distinct peaks in env
    hilb = hilbert_envelope(p_signal);
     
    if disp_ECG
        [e_signal, e_fs] = audioread(ECG_paths(i));
        disp(e_fs)
        e_dt = 1/e_fs; 
        e_t = 0:e_dt:(length(e_signal)*e_dt)-e_dt;

        M=e_fs;% define the window-length
        window=hanning(M); % define the window to use (Hanning)
        overlap=0.5; % define the overlap for windows
        [P,f]=spec2(e_signal,window,overlap,e_fs);% calculate the PSD
        figure();
        plot(f,P);% plot the result
        xticks(0:10:max(f));
        xlabel('frequency (Hz)');% label axes
        ylabel('PSD (microV^2/Hz)');
        p_title = strcat('PSD of ECG signal ', PCG_files(i).name)
        title(p_title);
       
        p_low={};
        p_low.N=2;
        p_low.type='butter';
        p_low.lphp='low';
        [b,a]=make_digital_filter(150, e_fs, p_low);
        e_signal = filtfilt(b, a, e_signal)
        
        p_high={};
        p_high.N=2;
        p_high.type='butter';
        p_high.lphp='high';
        [b,a]=make_digital_filter(1, e_fs, p_high);
        e_signal = filtfilt(b, a, e_signal);
        
        %e_signal = applyButterworthBandpassFilter(1, 20, 4, e_fs, e_signal);
        
        wo = 50/(e_fs/2);  
        bw = wo;
        [num,den] = iirnotch(wo,bw);
        e_signal = filtfilt(num, den, e_signal);
        
        windowWidth = 10; % Whatever you want.
        kernel = ones(windowWidth,1) / windowWidth;
        e_signal = filtfilt(kernel, 1, e_signal);
        
        
        e_signal = e_signal./max(abs(e_signal));
        e_signal = -e_signal; % Flip signal

        % display cleaned signal 
        figure();
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
        plot(t, p_signal, t, hilb);
    end

    if label_data
        good_rec = input(['Label this recording?'], 's')
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