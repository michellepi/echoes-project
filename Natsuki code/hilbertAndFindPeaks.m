function [heart_sounds,check] = hilbertAndFindPeaks(y,fs,window,ratio1,min_height,print,plot1,total_t)
% variable to check if working
check = 0;

% smoothing filter              
b = (1/window)*ones(1,window);  % parameter to make filter
a = 1;

%% hilbert
% applying hilbert energy envelope
hil = abs(hilbert(y));

% making an adaptive threshold
threshold = ratio1*max(hil);
l = length(hil);

% applying threshold
for j = 1:l
    if hil(j) > threshold(1,1)
        env_after_threshold(j) = hil(j);
    end
end
 
% applying smoothing to thresholded hilbert envelope
filtered = filter(b,a,env_after_threshold);

% normalising the hilbert data
filtered = filtered./max(filtered);

% find peaks
[~, Rt] = findpeaks(filtered, fs,'MinPeakHeight', min_height,'MinPeakDistance', 0.2);
Rt = Rt.';

% new time values
time = (1:length(filtered))./fs;

%% Plot result
% figure;
% plot(time,filtered)
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Plot hilbert transform');
% hold on

%% Find boundaries

% Average S1 and S2 duration in seconds
dur = 0.1;

for i = 1:length(Rt)
    bound_on(i) = Rt(i) - dur/2;
    bound_off(i) = Rt(i) + dur/2;
%     xline(Rt(i),'r','Linewidth', 0.5);
%     xline(bound_on(i),'g','Linewidth', 0.5);
%     xline(bound_off(i),'g','Linewidth', 0.5);
end

%% classify heart sounds as S1 or S2
[peaks_S1,peaks_S2,class] = S1_or_S2(Rt,print);

if print == 1
    for c = 1:length(class)
        if class(c) == 1
            fprintf('Peak %d: S1 \n', c);
        else
            fprintf('Peak %d: S2 \n', c);
        end
    end
end

%% Quit function and display error if bound_on or bound_off doesnt exist

% initialise heart_sounds matrix to return 0s if function quits
heart_sounds = zeros(length(class),3);

% if on time is before the start of time
for i = 1:length(bound_on)
    if bound_on(i) < 0
        bound_on(i) = 0.000001;
    else
        bound_on(i) = bound_on(i);
    end
end

% if off time is after the end of time
for i = 1:length(bound_off)
    if bound_off(i) > total_t
        bound_off(i) = total_t - 0.000001;
    else
        bound_off(i) = bound_off(i);
    end
end

% display error and quit if doesnt exist
try
    % reshaping
    bound_on = bound_on(bound_on>0);
    bound_off = bound_off(bound_off>0);
catch
    fprintf(2,'Cannot find necessary windows \n\n')
    return
end

%% creating matrix; 'class' which outputs '1' for S1 and '2' for S2
% creating 3 column matrix for each frame = 'heart_sounds' which is:
%      onset time     peak class     offset time

sz_on = size(bound_on); sz_off = size(bound_off);
if sz_on(1,1) == 1
    bound_on = transpose(bound_on);
end
if sz_off(1,1) == 1
    bound_off = transpose(bound_off);
end

% making sure everything is the right size
diff1_ = length(class) - length(bound_on);
diff2_ = length(class) - length(bound_off);
if length(class) > length(bound_on)
    bound_on = padarray(bound_on,diff1_,0,'post');
end
if length(class) > length(bound_off)
    bound_off = padarray(bound_off,diff2_,0,'post');
end

% make matrix
heart_sounds(:,1) = bound_on(:);
heart_sounds(:,2) = class(:);
heart_sounds(:,3) = bound_off(:);

%% Plot data with S1 and S2 visible
if plot1 == 1
    figure;
    range = 1:length(y);
    plot(range/fs,y)
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Plot hilbert transform with boundaries (S1 in Red)');
    hold on

    for i = 1:length(heart_sounds)
        if heart_sounds(i,2) == 1
            p(i) = (heart_sounds(i,3) - heart_sounds(i,1))/2 + heart_sounds(i,1);
            xline(heart_sounds(i,1),'r','Linewidth', 1.5)
            xline(heart_sounds(i,3),'r','Linewidth', 1.5)
        elseif heart_sounds(i,2) == 2
            p(i) = (heart_sounds(i,3) - heart_sounds(i,1))/2 + heart_sounds(i,1);
            xline(heart_sounds(i,1),'g','Linewidth', 1.5)
            xline(heart_sounds(i,3),'g','Linewidth', 1.5)
        end
    end
    xlim([0 7])
end


% checking variable
check = 1;

end