function mod_hs_peak = interpolatePeaks(first_classified_hs,target_interp_hs,start_or_end_seg,missingPeakIntv)
    
    if start_or_end_seg == 2 || start_or_end_seg == 1
        % find index of first classified heart sound right after end of missing peak interval
        first_classified_hs_idx1 = first_classified_hs(:)==missingPeakIntv(2);
        if isempty(find(first_classified_hs_idx1==true))
            mod_hs_peak=target_interp_hs;
            return
        end
        % find the first target peak we will use as reference for interpolation using peaks after missing peak interval 
        target_idx1 = find(round(target_interp_hs(:),4) >= missingPeakIntv(2)); 
        if ~isempty(target_idx1)
            target_idx1 = target_idx1(1);
        else 
            start_or_end_seg = 6;
        end
    end
    
    if start_or_end_seg == 2 || start_or_end_seg == 6
        % find index of first classified heart sound right before beginning of missing peak interval
        first_classified_hs_idx0 = find(first_classified_hs(:)==missingPeakIntv(1));
        if isempty(first_classified_hs_idx0)
            mod_hs_peak=target_interp_hs;
            return
        end
        % find the first target peak we will use as reference for interpolation using peaks before missing peak interval
        target_idx0 = find(round(target_interp_hs(:),4) <= missingPeakIntv(1)); 
        if ~isempty(target_idx0)
            target_idx0 = target_idx0(end);
        else 
            start_or_end_seg = 1;
        end
    end
    
    % if start_or_end_seg == 2, identify whether to inteprolate using peak before or after the inteval by nearest neighbour and changing its value to either 1 or 6
    if start_or_end_seg == 2
        diff1 = target_interp_hs(target_idx1) - missingPeakIntv(2);
        diff0 = missingPeakIntv(1) - target_interp_hs(target_idx0);

        if diff1 > diff0
            start_or_end_seg = 6;
        else
            start_or_end_seg = 1;
        end
    end
      
    if start_or_end_seg == 1 % interpolate using the heart sounds after the interval
        % calculate the difference between start of target and end of first classified hs
        diff_target_hs = target_interp_hs(target_idx1)-first_classified_hs(first_classified_hs_idx1);
        
        % create interpolated hs
        interp_hs = diff_target_hs + missingPeakIntv(1);
        
    elseif start_or_end_seg == 6 % inteprolate using the heart sounds before the interval    
        % calculate the difference between start of target and end of classified hs that doesn't contain the start of the missingPeakIntv
        diff_target_hs = target_interp_hs(target_idx0) - first_classified_hs(first_classified_hs_idx0-1);
        
        % create interpolated hs
        interp_hs = missingPeakIntv(1) + diff_target_hs;
        
    end
       
    % modify target hs array to include interpolated hs
    target_interp_hs(end+1) = interp_hs;
    mod_hs_peak = sortrows(target_interp_hs); % ensures heart sounds are in ascending order

end