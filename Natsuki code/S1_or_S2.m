function [peaks_S1,peaks_S2,class] = S1_or_S2(Rt1,print)

% pre-allocating matrices for speed
peaks_S1 = zeros(round(length(Rt1)),1);
peaks_S2 = zeros(round(length(Rt1)),1);
ind_S1   = zeros(round(length(Rt1)),1);
ind_S2   = zeros(round(length(Rt1)),1);

% calculate spaces_between matrix giving time differences between each peak
for c = 1:length(Rt1)-1
    space_b(c,1) = Rt1(c+1,1) - Rt1(c,1);
    spaces_between(c,1) = space_b(c,1);
end

% finds out if peaks are S1 or S2 [peaks excluding first and last peaks]
% ASSUMPTION: MINIMUM 3 PEAKS PER WINDOW
for c = 1:length(Rt1)
    if c < length(Rt1)-1
        if spaces_between(c,1) > spaces_between(c+1,1)
            ind_S1(c,1) = c+1;
            if ind_S1(c,1) ~= 0
                peaks_S1(c,1) = ind_S1(c,1);
            end
        else
            ind_S2(c,1) = c+1;
            if ind_S2(c,1) > 0
                peaks_S2(c,1) = ind_S2(c,1);
            end
        end
    end
end

% is the first peak S1 or S2
if ismember(2,ind_S1) == 1
    peaks_S2 = padarray(peaks_S2,1,1,'pre');
else
    peaks_S1 = padarray(peaks_S1,1,1,'pre');
end 

% is the last peak S1 or S2
if ismember(length(Rt1)-2,ind_S1) == 1
    peaks_S1(length(Rt1),1) = length(Rt1);
else
    peaks_S2(length(Rt1),1) = length(Rt1);
end

% deciding first peak
class = zeros(length(Rt1),1);
if ismember(1,peaks_S1) == 1
    class(1) = 1;
else
    class(1) = 2;
end

% deciding which peak is which taking into account that there might be missing peaks
for l = 2:length(Rt1)
    if ismember(l,peaks_S1) == 1
        class(l) = 1;
    elseif ismember(l,peaks_S1) == 0
        class(l) = 2;
    end
end

% printing out S1 and S2 values
if print == 1
    fprintf('\n');
    % print out classification results for S1
    for c = 1:length(class)
        if ismember(c,peaks_S1(c)) == 1
            fprintf('Peak %d is: S1 \n', c)
        end
    end

    fprintf('\n');
    % print out classification results for S2
    for d = 1:length(class)
        if ismember(d,peaks_S2(d)) == 1
            fprintf('Peak %d is: S2 \n', d)
        end
    end
end

end
