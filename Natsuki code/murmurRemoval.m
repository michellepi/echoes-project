function [y,y1,D,fs,total_t] = murmurRemoval(filename,low,high)

%% read in audio file and find ampltiude and time date
[y,fs] = audioread(filename);
y = y(:,1);
dt = 1/fs;
total_t = length(y)*dt;
t = (0:dt:(total_t)-dt)';

%% Filter 1: Bandpass
a = low;
b = high;
y1 = bandpass(y,[a b],fs,'Steepness',0.99);

% absolute values of bandpass filtered signal
y1 = abs(y1);
% normalised filtered signal
y1 = y1./max(y1);

%% Filter 2: Bandpass
a = 60;
b = 125;
y2 = bandpass(y,[a b],fs,'Steepness',0.9);

%% Thresholding coefficients

% threshold estimation
% 75th percentile of absolute values of signal, sorted in ascending order
val = abs(y2); val_sort = sort(val);
thresh = round(0.75*length(val_sort)); med75 = val_sort(thresh);

% calculate mean, median and variance of sorted absolute values
m = mean(val); MAD = median(val); v = var(val);

if med75 < v
    T = med75*(1-(v-med75));
elseif (med75 > v) && (med75 < m)
    T = med75;
elseif med75 > m
    T = med75 + (med75-m);
end

% threshold rescaling
% multiple scale dependent rescaling approach
noise_v = MAD/0.6745;
T = T*noise_v;

% applying threshold to coefficients
alpha = 1;
if med75 <= v
    beta = 1.3;
else
    beta = 1.4;
end
T1 = alpha*T; T2 = beta*T;

D = 0;
for i = 1:length(val)
    if val(i) > T2
        D(i) = val(i);
    elseif abs(val(i)) >= T1 && abs(val(i)) <= T2
        D(i) = (val(i)^3)/(T2^2);
    elseif val(i) < T1
        D(i) = 100;
    end
end
D = D(D>0); ind = find(D == 100); D(ind) = 0;

%% Scale larger amplitude data points
r = max(D) - min(D);
for p = 1:length(D)
    if D(p) > 0.5*r
        D(p) = D(p)/2;
    else
        D(p) = D(p);
    end
end

%% Normalising
D = D./max(abs(D));

%% Scale low amplitude noise data points
r = max(D) - min(D);
for p = 1:length(D)
    if D(p) < 0.1*r
        D(p) = 0;
    else
        D(p) = D(p);
    end
end

%% Plotting
% figure
% plot(t,y)
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Original Signal');
% 
% range1 = 1:length(D);
% figure;
% plot(range1./fs,D);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Processed Signal');

end

