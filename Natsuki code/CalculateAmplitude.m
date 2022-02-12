function [mean_of_amplitude12,mean_of_amplitude21,S1_to_S2_range_t,S2_to_S1_range_t,ratio] = CalculateAmplitude(heart_sounds,fs,orig,w)

%% checking if both heart sounds are visible
if ismember(2,heart_sounds(:,2)) == 0 || ismember(1,heart_sounds(:,2)) == 0
    S1_to_S2_range_t = zeros(20);
    S2_to_S1_range_t = zeros(20);
    mean_of_amplitude12 = zeros(20);
    mean_of_amplitude21 = zeros(20);
    ss = 0;
    return
end

%% Systolic range and amplitude calculations 

for e = 1:length(heart_sounds)-1
    if heart_sounds(e,2) == 1 && heart_sounds(e+1,2) == 2
        S1_to_S2_range_t(e,1) = heart_sounds(e,3);
        S1_to_S2_range_t(e,2) = heart_sounds(e+1,1);
    end
end

% turn the times into data point numbers
S1_to_S2_range = S1_to_S2_range_t(:,:).*fs;
siz = size(S1_to_S2_range);

% next loops can't deal with an array of size [1,1], pad with 0's
if siz(1) == 1
    S1_to_S2_range = padarray(S1_to_S2_range,1,0,'post');
end

% find amplitudes in S1 to S2 range
for f = 1:length(S1_to_S2_range)
    if S1_to_S2_range(f,1) > 0 && S1_to_S2_range(f,2) > 0 
        if heart_sounds(f,1) > 0
            a = round(S1_to_S2_range(f,1));
            b = round(S1_to_S2_range(f,2));

            amplitude = zeros(length(w),1);
            u = 1;
            for x = a+1:b
                amplitude(u) = w(x);
                u = u + 1;
            end
            amplitude = amplitude(amplitude>0);
            mean_of_amplitude12(f) = mean(amplitude);
        end
    end
end
mean_of_amplitude12 = mean_of_amplitude12(mean_of_amplitude12>0);

%% Diastolic range and amplitude calculations
S2_to_S1_range_t = zeros(length(heart_sounds),2);

for e = 1:length(heart_sounds)-1
    if heart_sounds(e,2) == 2 && heart_sounds(e+1,2) == 1
        S2_to_S1_range_t(e,1) = heart_sounds(e,3);
        S2_to_S1_range_t(e,2) = heart_sounds(e+1,1);
    end
end

% turn the times into data point numbers
S2_to_S1_range = S2_to_S1_range_t(:,:).*fs;
sz = size(S2_to_S1_range);

% find amplitudes in S2 to S1 range
for f = 1:sz(1)
    if S2_to_S1_range(f,1) > 0 && S2_to_S1_range(f,2) > 0 
        if heart_sounds(f,1) > 0
            a1 = round(S2_to_S1_range(f,1));
            b1 = round(S2_to_S1_range(f,2));

            amplitude2 = zeros(length(w),1);
            u = 1;
            for x = a1+1:b1
                amplitude2(u) = w(x);
                u = u + 1;
            end
        end
        amplitude2 = amplitude2(amplitude2>0);
        mean_of_amplitude21(f) = mean(amplitude2);
    end
end
mean_of_amplitude21 = mean_of_amplitude21(mean_of_amplitude21>0);

%% formatting of both S2-S1 and S1-S2

% replace negative values with 0s
for l = 1:length(mean_of_amplitude12)
    if mean_of_amplitude12(l) < 0 
        mean_of_amplitude12(l) = 0;
    end
end

for l = 1:length(mean_of_amplitude21)
    if mean_of_amplitude21(l) < 0 
        mean_of_amplitude21(l) = 0;
    end
end

while length(mean_of_amplitude12) > length(mean_of_amplitude21)
    mean_of_amplitude12(end) = [];
end
while length(mean_of_amplitude21) > length(mean_of_amplitude12)
    mean_of_amplitude21(end) = [];
end

%% Calculating individual cycle S/D ratios
ratio = zeros(100,1);

% calculate each systolic to diastolic ratio
for i = 1:length(mean_of_amplitude12)
    ratio(i) = mean_of_amplitude12(i)/mean_of_amplitude21(i);
end

% make sure ratio is column-orientated
s = size(ratio);
if s(2)> 1
    ratio = ratio';
end

end