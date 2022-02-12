function [ind_S1, ind_S2, missingPeakIntv, missingPeakIntvWidth, rm_idx] = classifyPeaks(joined_pks_loc,joined_ref_idx,sysTimeIntv,heartRate)

% pre-allocating matrices for speed
peaks_S1 = zeros(round(length(joined_pks_loc)),1);
peaks_S2 = zeros(round(length(joined_pks_loc)),1);
ind_S1   = zeros(round(length(joined_pks_loc)),1);
ind_S2   = zeros(round(length(joined_pks_loc)),1);
missingPeakIntv = []; missingPeakIntvWidth = []; prev_peak = []; rm_idx = [];

diaTimeIntv = (60/heartRate-sysTimeIntv);

% calculate difference between peaks
intv = diff(joined_pks_loc); 

% identify peaks that are too close to one another (one of them is an artefact)
small_intv = find(intv<0.85*sysTimeIntv);
if ~isempty(small_intv)
    % remove appropriate peaks and re-calculate new intv
    for i=1:length(small_intv)
        if small_intv(i) > 1
            if intv(small_intv(i)+1) < intv(small_intv(i)-1) && (joined_pks_loc(small_intv(i)+1)-joined_pks_loc(small_intv(i)-1)) < 1.25*diaTimeIntv
                rm_idx(end+1) = small_intv(i);
            elseif (joined_pks_loc(small_intv(i)+2)-joined_pks_loc(small_intv(i))) < 1.25*diaTimeIntv
                rm_idx(end+1) = small_intv(i)+1;
            end
        elseif small_intv(i) == 1 && (joined_ref_idx(2,2)-joined_ref_idx(1,1) < sysTimeIntv)
            rm_idx(end+1) = small_intv(i);
        end
    end
end
joined_pks_loc(rm_idx) = []; joined_ref_idx(rm_idx,:) = [];
intv = diff(joined_pks_loc);

rm_pk = find(intv>1.25*diaTimeIntv);
while ~isempty(rm_pk) && rm_pk(1) ==1
    joined_pks_loc(rm_pk(1)) = []; joined_ref_idx(rm_pk(1),:) = [];
    intv = diff(joined_pks_loc);
    rm_pk = find(intv>1.25*diaTimeIntv);
end

% identify missing peaks
missingPeakIntv_idx=find(intv > 1.25*diaTimeIntv);
if ~isempty(missingPeakIntv_idx)
    for i=1:length(missingPeakIntv_idx)
        missingPeakIntv(i,:) = [joined_pks_loc(missingPeakIntv_idx(i)), joined_pks_loc(missingPeakIntv_idx(i)+1)];
        missingPeakIntvWidth(i,:) = [joined_ref_idx(missingPeakIntv_idx(i),2), joined_ref_idx(missingPeakIntv_idx(i)+1,1)];
    end
end

c=1;
while c < size(intv,1)
    if intv(c,1) > intv(c+1,1) 
        if intv(c+1) > 1.25*(60/heartRate-sysTimeIntv) && intv(c,1) > 1.25*diaTimeIntv
            % if this case is true, there is a missing peak
            if prev_peak(end) == true
                ind_S1(c,1) = c+1; prev_peak(end+1) = true;
            else
                ind_S2(c,1) = c+1; prev_peak(end+1) = false;
            end
        elseif intv(c,1) > 1.25*(60/heartRate-sysTimeIntv) && intv(c+1,1) < 1.25*diaTimeIntv
            if prev_peak(end) == true
                ind_S1(c,1) = c+1; prev_peak(end+1) = true;
            else
                ind_S2(c,1) = c+1; prev_peak(end+1) = false;
            end
        else
            ind_S1(c,1) = c+1; prev_peak(end+1) = true;
        end
    elseif intv(c,1) < intv(c+1,1) 
        if intv(c+1) > 1.25*(60/heartRate-sysTimeIntv) && intv(c,1) > 1.25*diaTimeIntv
            if prev_peak(end) == true
                ind_S1(c,1) = c+1; prev_peak(end+1) = true;
            else
                ind_S2(c,1) = c+1; prev_peak(end+1) = false;
            end
        elseif intv(c+1) > 1.25*(60/heartRate-sysTimeIntv) && intv(c,1) < 1.25*diaTimeIntv
            if prev_peak(end) == true
                ind_S2(c,1) = c+1; prev_peak(end+1) = false;
            else
                ind_S1(c,1) = c+1; prev_peak(end+1) = true;
            end
        else 
            ind_S2(c,1) = c+1; prev_peak(end+1) = false;
        end
    end
    c = c+1;
end

% classify first peak
if ismember(2,ind_S1)
    if intv(1) > 1.25*diaTimeIntv
        if prev_peak(1) == true
            ind_S1 = [1; ind_S1];
        else
            ind_S2 = [1; ind_S2];
        end
    else
        ind_S2 = [1; ind_S2];
    end
elseif ismember(2,ind_S2)
    if intv(1) > 1.25*diaTimeIntv
        if prev_peak(1) == true
            ind_S1 = [1; ind_S1];
        else
            ind_S2 = [1; ind_S2];
        end
    else
       ind_S1 = [1; ind_S1];
    end
end

% classify last peak
if ismember(length(intv),ind_S1)
    if intv(end) > 1.25*diaTimeIntv
        if prev_peak(end) == true
            ind_S1(end+1) = length(intv)+1;
        else
            ind_S2(end+1) = length(intv)+1;
        end
    else
        ind_S2(end+1) = length(intv)+1;
    end
elseif ismember(length(intv),ind_S2)
    if intv(end) > 1.25*diaTimeIntv
        if prev_peak(end) == true
            ind_S1(end+1) = length(intv)+1;
        else
            ind_S2(end+1) = length(intv)+1;
        end
    else
       ind_S1(end+1) = length(intv)+1;
    end
end

end