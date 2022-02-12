function [s1_refined, s2_refined, s1_peaks, s2_peaks, missingHS] = finds1s2_4(segmentedPeaks,fs,t,heartRate,sysTimeIntv)
%% testing function
% [segmentedPeaks, fs, t, audio_resampled, t2, sysTimeIntv, heartRate] = pre_processAudio2(1,1);

%% figures
% postFindPeaks = figure('Name', 'Peaks identified by FindPeaks');
% postPeakRefine = figure('Name', 'After peak refinement');

%% global variables
dt = 1/fs;                                                                  % time length of one sample
L = size(t,2);                                                              % number of segments
ref_idx = cell(L,1);                                                        % cells to store peak indices identified by findPeaks
pks_stored = cell(L,1);                                                     % cell stores peak values 
pks_loc_stored = cell(L,1);                                                 % cell stores peak locations
joined_ref_idx = [];                                                        % stores joined peak indices 
joined_pks = [];                                                            % stores joined peak values
joined_pks_loc = [];                                                        % stores joined peak locations
missingPeakIntv = zeros(50,2);                                              % intervals indicating missing peaks
pad_zeros = zeros(50,1);                                                    % number of zeros to pad array with so findpeaks works well
intv_thresh = 0.1;                                                          % peaks of heart sounds cannot be closer than this threshold (in seconds)
missingHS = zeros(1,2);                                                     % array holds number of missing peaks of s1 and s2 sounds

%% begin s1 and s2 classification by segment
for z=1:L
    %% first refine peaks using systolic time interval and heart rate to remove unwanted peaks
    % find current peak locations
    [idx0, idx1] = findPeakLoc(segmentedPeaks(:,z));
    
    % to compensate for matlab indices starting at 1
    idx0 = idx0 - 1; idx1 = idx1 - 1;
    
    % find peaks that are separated by at least the systolic time interval 
    padded_segment = [pad_zeros; segmentedPeaks(:,z); pad_zeros];
    padded_t_array = t(1,z):dt:t(1,z)+(length(padded_segment)*dt)-dt; 
    [pks, pks_loc] = findpeaks(padded_segment, padded_t_array,'MinPeakDistance',sysTimeIntv*0.8);
    
    % re-adjust peak locations so they take into account padding
    pks_loc = pks_loc - length(pad_zeros)*dt;

    %%  displaying identified by findPeaks-----------------------------------------------------------
%     figure(postFindPeaks);
%     subplot(2,3,z);
%     plot(t(:,z), segmentedPeaks(:,z));hold on;
%     plot(pks_loc,pks,'r*');
%     xlabel('Time(s)');ylabel('Amplitude'); 
%     title('FindPeaks');
%  
    %% only store peak locations that have been identified --------------------------------------
    % intialise array to store identified peak intervals
    ref_idx0 = zeros(length(pks_loc),2);
    
    % storing refined start and end time values of peaks
    for i = 1:length(pks_loc)
        for j = 1:length(idx0)
            if (round(pks_loc(i),4) >= round((idx0(j)*dt)+(z-1),4)) && (round(pks_loc(i),4) <= round((idx1(j)*dt)+(z-1),4))
                ref_idx0(i,1) = (idx0(j)*dt) + (z-1); 
                ref_idx0(i,2) = (idx1(j)*dt) + (z-1); 
                break
            end
        end
    end
    
    % store identified peak indices so they can be joined together
    ref_idx{z,1} = ref_idx0;
    pks_stored{z,1} = pks;
    pks_loc_stored{z,1} = pks_loc';
end

%% join all segments together and remove overlapping peaks
for z = 1:L-1
    rm_idx = []; % stores peak indices that should be removed
    ref_idx_1 = cell2mat(ref_idx(z)); ref_idx_2 = cell2mat(ref_idx(z+1));
    ref_idx_cat = cat(1,ref_idx_1,ref_idx_2); 
    ref_idx_cat_sorted = sortrows(ref_idx_cat);
    
    % do the same for peak values and peak loc cell arrays
    pks_1 = cell2mat(pks_stored(z)); pks_2 = cell2mat(pks_stored(z+1));
    pks_cat = cat(1,pks_1,pks_2);
    pks_loc1 = cell2mat(pks_loc_stored(z)); pks_loc2 = cell2mat(pks_loc_stored(z+1));
    pks_loc_cat = cat(1,pks_loc1,pks_loc2);
    
    % calculate mean peak width
    mean_pk_width = mean(ref_idx_cat(2,:)-ref_idx_cat(1,:));
    
    % check for overlapping peaks in segments
    for i = 1:size(ref_idx_cat,1)-1
        [remove_peak, overlap_TF] = correctForOverlappingPeaks(ref_idx_cat_sorted(i,:), ref_idx_cat_sorted(i+1,:),mean_pk_width);
        if overlap_TF == false
            continue;
        else
            rm_idx(end+1,:) = remove_peak;
        end
    end
    for i = 1:size(rm_idx,1)
        pks_cat(find(ref_idx_cat(:,1) == rm_idx(i,1))) = [];
        pks_loc_cat(find(ref_idx_cat(:,1) == rm_idx(i,1))) = [];
        ref_idx_cat(find(ref_idx_cat(:,1) == rm_idx(i,1)),:) = [];
    end
    joined_ref_idx = cat(1,joined_ref_idx,ref_idx_cat);
    joined_pks = cat(1,joined_pks,pks_cat);
    joined_pks_loc = cat(1,joined_pks_loc,pks_loc_cat);
    
end
% remove duplicates in joined segments
[joined_ref_idx unique_idx] = unique(round(joined_ref_idx,4),'rows');
joined_pks = joined_pks(unique_idx);
joined_pks_loc = joined_pks_loc(unique_idx);

%% deleting overlapping peaks after joining segments together
rm_idx = []; 
for i = 1:size(joined_ref_idx,1)-1
    % calculate mean peak width
    mean_pk_width = mean(joined_ref_idx(:,2) - joined_ref_idx(:,1));
    
    [remove_peak, overlap_TF] = correctForOverlappingPeaks(joined_ref_idx(i,:), joined_ref_idx(i+1,:), mean_pk_width);
    if overlap_TF == false
        continue;
    else
        rm_idx(end+1,:) = remove_peak;
    end
end
for j = 1:size(rm_idx,1) 
    joined_pks(ismember(joined_ref_idx,rm_idx(j,:),'rows'))=[];
    joined_pks_loc(ismember(joined_ref_idx,rm_idx(j,:),'rows'))=[];
    joined_ref_idx(ismember(joined_ref_idx,rm_idx(j,:),'rows'),:)=[];
end
% % remove duplicates in joined segments
[joined_ref_idx unique_idx] = unique(round(joined_ref_idx,4),'rows');
joined_pks = joined_pks(unique_idx);
joined_pks_loc = joined_pks_loc(unique_idx);

%% removal of peaks right next to each other (ie <= 100ms)
[joined_ref_idx, joined_pks, joined_pks_loc] = removeClosePeaks(joined_ref_idx, joined_pks, joined_pks_loc, intv_thresh);

% figure(15); 
% plot(t2,audio_resampled,'k'); hold on; plot(joined_pks_loc,joined_pks,'*r'); hold on;
% for i=1:size(joined_ref_idx,1)
%     hold on; xline(joined_ref_idx(i,1),'r');
%     hold on; xline(joined_ref_idx(i,2),'r');
% end

%% classifying peaks
[s1_refined, s2_refined, missingPeakIntv, missingPeakIntvWidth, rm_idx] = classifyPeaks(joined_pks_loc,joined_ref_idx,sysTimeIntv,heartRate);

% updating identified peaks if any have to be removed
if ~isempty(rm_idx)
    joined_pks_loc(rm_idx) = []; joined_ref_idx(rm_idx,:) = []; joined_pks(rm_idx) = [];
end
s1_idx = find(s1_refined~=0); s2_idx = find(s2_refined~=0);

% these contain the locations of the heart sound peaks
s1_peaks = joined_pks_loc(s1_refined(s1_idx)); s2_peaks = joined_pks_loc(s2_refined(s2_idx));
% these contain heart sound time intervals
s1_refined = joined_ref_idx(s1_refined(s1_idx),:); s2_refined = joined_ref_idx(s2_refined(s2_idx),:);

% figure(21); 
% plot(t,segmentedPeaks,'LineWidth',1); hold on; plot(joined_pks_loc,joined_pks,'*r'); hold on;
% plot(t2,audio_resampled,'k'); hold on;
% 
% s1_idx = find(s1_refined~=0); s2_idx = find(s2_refined~=0);
% for i = 1:size(s1_idx,1)
%     hold on; xline(joined_ref_idx(s1_refined(s1_idx(i)),1),'r');
%     hold on; xline(joined_ref_idx(s1_refined(s1_idx(i)),2),'r');
% end
% for i = 1:size(find(s2_refined~=0),1)
%     hold on; xline(joined_ref_idx(s2_refined(s2_idx(i)),1),'g');
%     hold on; xline(joined_ref_idx(s2_refined(s2_idx(i)),2),'g');
% end

%% nearest-neighbour interpolation of missing peaks 
% removing duplicate intervals
missingPeakIntv = unique(missingPeakIntv,'rows');
missingPeakIntvWidth = unique(missingPeakIntvWidth,'rows');

if ~isempty(missingPeakIntvWidth)
    for i = 1:size(missingPeakIntv,1)
        % identify classified peak after interval
        first_classified_hs = find(s1_refined==missingPeakIntvWidth(i,2));

        if ~isempty(first_classified_hs)
            % this means peak directly after it is an s1 peak therefore we interpolate an s2 peak into the interval
            mod_hs = interpolatePeaksWidth(s1_refined,s2_refined,2,missingPeakIntvWidth(i,:));
            mod_hs_peak = interpolatePeaks(s1_peaks,s2_peaks,2,missingPeakIntv(i,:));
            % replace modified heart sound into s2
            s2_refined = mod_hs; s2_peaks = mod_hs_peak;
            % identify missing heart sound as s2
            missingHS(2) = missingHS(2)+1;
        else
            mod_hs = interpolatePeaksWidth(s2_refined,s1_refined,2,missingPeakIntvWidth(i,:));
            mod_hs_peak = interpolatePeaks(s2_peaks,s1_peaks,2,missingPeakIntv(i,:));
            % replace modified heart sound into s1
            s1_refined = mod_hs; s1_peaks = mod_hs_peak;
            % identify missing heart sound as s1
            missingHS(1) = missingHS(1)+1;
        end
    end
    
end
% remove duplicates in refined s1 and s2
s1_refined = unique(round(s1_refined,4),'rows'); s2_refined = unique(round(s2_refined,4),'rows');
s1_peaks = unique(round(s1_peaks,4)); s2_peaks = unique(round(s2_peaks,4));

%% displaying final results
% figure(21); 
% plot(t,segmentedPeaks,'LineWidth',1); hold on; 
% plot(t2,audio_resampled,'k'); hold on;
% 
% for i = 1:size(s1_refined,1)
%     hold on; xline(s1_refined(i,1),'r');
%     hold on; xline(s1_peaks(i),'-m');
%     hold on; xline(s1_refined(i,2),'r');
% end
% for i = 1:size(s2_refined,1)
%     hold on; xline(s2_refined(i,1),'g');
%     hold on; xline(s2_peaks(i),'-b');
%     hold on; xline(s2_refined(i,2),'g');
% end
  
end






%% START OF FUNCTIONS -----------------------------------------------------
%% ------------------------------------------------------------------------
function [peaks, peak_vals, peak_locs] = removeClosePeaks(peaks,peak_vals,peak_locs,intv_thresh)
    rm_peak = [];
    
    % calculate differences between peaks
    intv_of_pks = peaks(2:end,1) - peaks(1:end-1,2);
    pk_idx = find(abs(intv_of_pks) < intv_thresh);

    for i = 1:length(pk_idx)
        %compare peak widths and choose to remove narrower peak
        width1 = peaks(pk_idx(i),2) - peaks(pk_idx(i),1);
        width2 = peaks(pk_idx(i)+1,2) - peaks(pk_idx(i)+1,1);
        if width1 > width2
            rm_peak(end+1) = pk_idx(i)+1;
        else
            rm_peak(end+1) = pk_idx(i);
        end 
    end
    % remove peaks
    peaks(rm_peak,:) = [];
    peak_vals(rm_peak,:) = [];
    peak_locs(rm_peak,:) = [];
end

%% ------------------------------------------------------------------------
function [remove_peak, overlap_TF] = correctForOverlappingPeaks(peak1, peak2,mean_pk_width)
% function corrects for overlapping peaks after segments have been joined together. This works on only 2 consecutive peaks.
% outputs:
%   - remove_peak: the value of the peak that should be removed, not the index of it
%   - overlap_TF: booloean values which indicates is the two peaks overlapped

    overlap_TF = false; remove_peak = 0;
    if (peak1(2) > peak2(1) && peak1(1) <= peak2(1)) || (peak2(2) > peak1(1) && peak2(1) <= peak1(1)) || (peak1(1) < peak2(1) && peak1(2) > peak2(2)) || (peak2(1) < peak1(1) && peak2(2) > peak1(2))
        % calculate widths and return peak with smaller width to be removed
        width1 = peak1(2) - peak1(1);
        width2 = peak2(2) - peak2(1);
        if width1 > width2 || width1 == width2
            if width1 > mean_pk_width
                remove_peak = peak1;
            else
                remove_peak = peak2;
            end
        elseif width2 > width1
            if width2 > mean_pk_width
                remove_peak = peak2;
            else
                remove_peak = peak1;    
            end
        end
        overlap_TF = true;
    end
end

%% ------------------------------------------------------------------------
function [mod_hs] = interpolatePeaksWidth(first_classified_hs,target_interp_hs,start_or_end_seg,missingPeakIntv)
    % inputs:
    %   - first_classified_hs: the classified heart sound array belonging to the first heart sound of the missing peak interval
    %   - target_interp_hs: the heart sound array belonginging to the second heart sounds identified after the missing interval
    %   - start_or_end_seg: integer value of either 1, 6 or 2 where 1 and 6 identifies whether it is the first or last segment and 2 identifies that it is neither and so can use either last identified HS (6) or the heart sound after it (1)
    
    if start_or_end_seg == 2 || start_or_end_seg == 1
        % find index of first classified heart sound right after end of missing peak interval
        first_classified_hs_idx1 = first_classified_hs(:,1)==missingPeakIntv(2);
        if isempty(find(first_classified_hs_idx1==true))
            mod_hs=target_interp_hs;
            return
        end
        % find the first target peak we will use as reference for interpolation using peaks after missing peak interval 
        target_idx1 = find(round(target_interp_hs(:,1),4) >= missingPeakIntv(2)); 
        if ~isempty(target_idx1)
            target_idx1 = target_idx1(1);
        else 
            start_or_end_seg = 6;
        end
    end
    
    if start_or_end_seg == 2 || start_or_end_seg == 6
        % find index of first classified heart sound right before beginning of missing peak interval
        first_classified_hs_idx0 = find(first_classified_hs(:,2)==missingPeakIntv(1));
        if isempty(first_classified_hs_idx0)
            mod_hs=target_interp_hs;
            return
        end
        % find the first target peak we will use as reference for interpolation using peaks before missing peak interval
        target_idx0 = find(round(target_interp_hs(:,2),4) <= missingPeakIntv(1)); 
        if ~isempty(target_idx0)
            target_idx0 = target_idx0(end);
        else 
            start_or_end_seg = 1;
        end
    end
    
    % if start_or_end_seg == 2, identify whether to inteprolate using peak before or after the inteval by nearest neighbour and changing its value to either 1 or 6
    if start_or_end_seg == 2
        diff1 = target_interp_hs(target_idx1,1) - missingPeakIntv(2);
        diff0 = missingPeakIntv(1) - target_interp_hs(target_idx0,2);

        if diff1 > diff0
            start_or_end_seg = 6;
        else
            start_or_end_seg = 1;
        end
    end
      
    if start_or_end_seg == 1 % interpolate using the heart sounds after the interval
        % find the width of the target hs for interpolation
        width_target_hs = target_interp_hs(target_idx1,2)-target_interp_hs(target_idx1,1);
        
        % calculate the difference between start of target and end of first classified hs
        diff_target_hs = target_interp_hs(target_idx1,1)-first_classified_hs(first_classified_hs_idx1,2);
        
        % create interpolated hs
        interp_hs = [diff_target_hs, diff_target_hs+width_target_hs] + missingPeakIntv(1);
        
    elseif start_or_end_seg == 6 % inteprolate using the heart sounds before the interval    
        % find the width of the target hs for interpolation
        width_target_hs = target_interp_hs(target_idx0,2)-target_interp_hs(target_idx0,1);
        
        % calculate the difference between start of target and end of classified hs that doesn't contain the start of the missingPeakIntv
        diff_target_hs = target_interp_hs(target_idx0,1) - first_classified_hs(first_classified_hs_idx0-1,2);
        
        % create interpolated hs
        interp_hs = missingPeakIntv(1) + [diff_target_hs+width_target_hs diff_target_hs];
        
    end
       
    % modify target hs array to include interpolated hs
    target_interp_hs(end+1,:) = interp_hs;
    mod_hs = sortrows(target_interp_hs); % ensures heart sounds are in ascending order

end

