function [noise_level noise_score] = calculateNoiseLevel(signal)
% output: 
%   - noise_level: calculated noise level of recording
%   - noise_score: quantified noise level of recording

%% threshold estimation
% 75th percentile of absolute values of signal, sorted in ascending order
val = abs(signal); val_sort = sort(val);
thresh = round(0.75*length(val_sort)); med75 = val_sort(thresh);

% calculate mean, median and variance of sorted absolute values
m = mean(val);
MAD = median(val);
v = var(val);

if med75 < v
    T = med75*(1-(v-med75));
    noise_level = 'low';
elseif (med75 > v) && (med75 < m)
    T = med75;
    noise_level = 'medium';
elseif med75 > m
    T = med75 + (med75-m);
    noise_level = 'high';
end

%% quantifying noise level
diff_v = (med75 - v)*100/med75;
diff_m = (med75 - m)*100/med75;
noise_score = diff_v + diff_m;
disp(['Noise score: ',num2str(noise_score)]);

end