function [mean_t,number,s1_m,s1_std,s2_m,s2_std] = AcceptableWindows(heart_sounds)
% tolerance
tol = 0.15;

%% initialise variables
S1_start_t = zeros(length(heart_sounds));
S1_end_t = zeros(length(heart_sounds));
S2_start_t = zeros(length(heart_sounds));
S2_end_t = zeros(length(heart_sounds));

%% calculations of times
% create matrix of S1 times
for n = 1:length(heart_sounds)
    if heart_sounds(n,2) == 1
        S1_start_t(n) = heart_sounds(n,1); S1_end_t(n) = heart_sounds(n,3);
    end
end
S1_start_t = S1_start_t(S1_start_t>0)'; S1_end_t = S1_end_t(S1_end_t>0)';

% calculate the time between S1s
for n = 2:length(S1_start_t)
    S1_s_diff(n) = S1_start_t(n) - S1_start_t(n-1);
end
for n = 2:length(S1_end_t)
    S1_e_diff(n) = S1_end_t(n) - S1_end_t(n-1);
end
S1_s_diff = S1_s_diff(S1_s_diff>0.5)'; S1_e_diff = S1_e_diff(S1_e_diff>0.5)';
S1_s_diff = S1_s_diff(S1_s_diff<1.2); S1_e_diff = S1_e_diff(S1_e_diff<1.2);

% create matrix of S2 times
for n = 1:length(heart_sounds)
    if heart_sounds(n,2) == 2
        S2_start_t(n) = heart_sounds(n,1); S2_end_t(n) = heart_sounds(n,3);
    end
end
S2_start_t = S2_start_t(S2_start_t>0)'; S2_end_t = S2_end_t(S2_end_t>0)';

% calculate the time between S2s
for n = 2:length(S2_start_t)
    S2_s_diff(n) = S2_start_t(n) - S2_start_t(n-1);
end
for n = 2:length(S2_end_t)
    S2_e_diff(n) = S2_end_t(n) - S2_end_t(n-1);
end
S2_s_diff = S2_s_diff(S2_s_diff>0.5)';S2_e_diff = S2_e_diff(S2_e_diff>0.5)';
S2_s_diff = S2_s_diff(S2_s_diff<1.2); S2_e_diff = S2_e_diff(S2_e_diff<1.2);

%% Look for outliers with tolerance of +-15%
% find the most common rounded time difference of each metric
rounded1 = round(S1_s_diff,1); mc(1) = mode(rounded1,'all');
rounded2 = round(S1_e_diff,1); mc(2) = mode(rounded2,'all');
rounded3 = round(S2_s_diff,1); mc(3) = mode(rounded3,'all');
rounded4 = round(S2_e_diff,1); mc(4) = mode(rounded4,'all');

% check if all boundary differences were found
if isempty(S1_s_diff) == 1 || isempty(S1_e_diff) == 1 || isempty(S2_s_diff) == 1 || isempty(S2_e_diff) == 1
    number = 0;
    mean_t = 0;
    return
end

% S1 - S1 start time differences
for l = 1:length(S1_s_diff)
    if S1_s_diff(l) > mc(1)-tol*mc(1) && S1_s_diff(l) < mc(1)+tol*mc(1)
        S1_s_all(l) = S1_s_diff(l);
    end
end
S1_s = S1_s_all(S1_s_all>0); count(1) = length(S1_s); m(1) = mean(S1_s); v(1) = var(S1_s);

% S1 - S1 end time differences
for l = 1:length(S1_e_diff)
    if S1_e_diff(l) > mc(2)-tol*mc(2) && S1_e_diff(l) < mc(2)+tol*mc(2) 
        S1_e_all(l) = S1_e_diff(l);
    end
end
S1_e = S1_e_all(S1_e_all>0); count(2) = length(S1_e); m(2) = mean(S1_e); v(2) = var(S1_e);

% S2 - S2 start time differences
for l = 1:length(S2_s_diff)
    if S2_s_diff(l) > mc(3)-tol*mc(3) && S2_s_diff(l) < mc(3)+tol*mc(3)
        S2_s_all(l) = S2_s_diff(l);
    end
end
S2_s = S2_s_all(S2_s_all>0); count(3) = length(S2_s); m(3) = mean(S2_s); v(3) = var(S2_s);

% S2 - S2 end time differences
for l = 1:length(S2_e_diff)
    if S2_e_diff(l) > mc(4)-tol*mc(4) && S2_e_diff(l) < mc(4)+tol*mc(4)
        S2_e_all(l) = S2_e_diff(l);
    end
end
S2_e = S2_e_all(S2_e_all>0); count(4) = length(S2_e); m(4) = mean(S2_e); v(4) = var(S2_e);

%% metric calculations
% mean length of all 'acceptable windows'
mean_t = mean(m(1:4));

% mean and std of S1 intervals
s1_m = mean(m(1:2));
s1_var = mean(v(1:2)); s1_std = sqrt(s1_var);

% mean and std of S2 intervals
s2_m = mean(m(3:4));
s2_var = mean(v(3:4)); s2_std = sqrt(s2_var);

% number of acceptable windows
number = floor(mean(count(:)));

end
