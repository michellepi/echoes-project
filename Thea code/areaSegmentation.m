function segmentedPeaks = areaSegmentation(thresh_hilbEnergy, fs, t)
% area segmentation process based on paper by Sharma et al (referring to Sharma et al 
%% for testing:
% [thresh_hilbEnergy, t, fs] = pre_processAudio2(3,3);

%% figures 
%    areaSeg = figure('Name','After Area Segmentation');
%    peakSep = figure('Name','After Peak Separation');
%    
%% global variables
    L = size(t,2);                                                          % the number of segments being used for the recording
    area = cell(L,1);                                                       % cell array to store area values for peaks fo segments
    segmentedPeaks = thresh_hilbEnergy;                                     % segmentedPeaks array to be manipulated for output
    idx = cell(L,2);                                                        % stores indices of peak locations
    dt = 1/fs;                                                              % the time interval length of a single sample
    theta = zeros(L,1);                                                     % array stores the theta values of each segment
    theta(1) = 0.05;                                                        % initialise first theta value to be used for first segment
    intv_thresh = 0.1/dt;                                                   % the interval in number of samples rather than time

%% for each segment: -----------------------------------------------------------------------
    for z = 1:L
        % find locations of peaks
        [idx0 idx1 thresh_hilbEnergy(:,z)] = findPeakLoc(thresh_hilbEnergy(:,z));
        idx{z,1} = idx0; idx{z,2} = idx1;
        
        % calculate area of each peak
        tmp_area = calcPeakArea(idx0,idx1,thresh_hilbEnergy(:,z),t(:,z));
        area{z,1} = tmp_area;
        
        % unless using first segment, change theta according to previous segment value
        if z ~=1
            % change theta according to mean of 3 peaks with highest area in previous segment
            [max_areas max_areas_idx] = maxk(cell2mat(area(z-1)),3);
            idx0_old = cell2mat(idx(z-1,1)); idx1_old = cell2mat(idx(z-1,2));
            theta(z) = mean(idx1_old(max_areas_idx)-idx0_old(max_areas_idx))*dt;
        end
        
        % array of peak widths
        dx = (idx1-idx0).*dt;
        % separating peaks according to criteria
        for i = 1:length(dx)
            split_idx=[];
            if (dx(i) > 0) && (dx(i) <= (2.5*theta(z)))
                continue; % keep peak the same
            elseif (dx(i) > (2.5*theta(z))) && (dx(i) <= (3*theta(z)))
                % split peak into 2 parts
                split_idx = splitPeakIdx(idx0(i),idx1(i),2);     
            elseif (dx(i) > (3*theta(z))) && (dx(i) <= (4*theta(z)))
                % split peak into 3 parts
                split_idx = splitPeakIdx(idx0(i),idx1(i),3);       
            elseif (dx(i) > (4*theta(z))) && (dx(i) <= (5*theta(z)))
                % split peak into 4 parts
                split_idx = splitPeakIdx(idx0(i),idx1(i),4);
            elseif (dx(i) > (5*theta(z))) && (dx(i) <= (6*theta(z)))
                % split peak into 5 parts
                split_idx = splitPeakIdx(idx0(i),idx1(i),5);
            elseif (dx(i) > (6*theta(z))) && (dx(i) <= (7*theta(z)))
                % split peak into 5 parts
                split_idx = splitPeakIdx(idx0(i),idx1(i),6);
            end
            segmentedPeaks(split_idx,z) = 0;
        end
 
        % update peak locations and area measurements
        [idx0 idx1 segmentedPeaks(:,z)] = findPeakLoc(segmentedPeaks(:,z));
        idx{z,1} = idx0; idx{z,2} = idx1;                                                
        
        tmp_area = calcPeakArea(idx0,idx1,thresh_hilbEnergy(:,z),t(:,z));
        area{z,1} = tmp_area;
        
        
        %% displaying separated peaks ------------------------------------------------------------------
%         figure(peakSep);
%         subplot(1,L,z);
%         plot(t(:,z), segmentedPeaks(:,z));
%         xlabel('Time');ylabel('Amplitude');
%         title(['Segment:',num2str(z),' After Peak Separation']);

        %% removing peaks that are within 100ms to each other ------------------------------------------
        intv_width_TF = calcIntvWidthBelowThresh(idx0,idx1,intv_thresh,segmentedPeaks(:,z),t(:,z));
        while ~isempty(intv_width_TF)
            
            % delete peak corresponding to first element of intv_width_TF
            segmentedPeaks(idx0(intv_width_TF):idx1(intv_width_TF),z) = 0;
            
            % find new peak locations
            [idx0 idx1 segmentedPeaks(:,z)] = findPeakLoc(segmentedPeaks(:,z));
            
            % find new intv_width_TF
            intv_width_TF = calcIntvWidthBelowThresh(idx0,idx1,intv_thresh,segmentedPeaks(:,z),t(:,z));
            
        end
        
        %% -----------------------------------------------------------------------------------------------
        % displaying peaks after full area segmentation
%         figure(areaSeg);
%         subplot(1,L,z);
%         plot(t(:,z), segmentedPeaks(:,z));
%         xlabel('Time');ylabel('Amplitude'); 
%         title(['Segment:',num2str(z),' After Area Segmentation']);


    end   
end


%% ----------------------------------------------------------------------------------------

function tmp_area = calcPeakArea(idx0,idx1,hilbEnergy,t)
% function calculates the area under a peak
    tmp_area = zeros(length(idx0),1);
    for i = 1:length(idx0)
        tmp_area(i) = trapz(t(idx0(i):idx1(i)),hilbEnergy(idx0(i):idx1(i)));
    end
end 


%% ----------------------------------------------------------------------------------------

function split_idx = splitPeakIdx(idx0,idx1,num_split)
% function returns indices in which to split peaks by
    
    if num_split == 2
        % split peaks from the sides rather than even numbers
        cut_idx = (idx1-idx0)*0.2;
        split_idx = [idx0+floor(cut_idx/2) idx1-floor(cut_idx/2)];
        return
    else
        split_idx = round(linspace(idx0,idx1,num_split+1));
        split_idx(1) = [];
        split_idx(end) = [];
    end
    
end



%% ----------------------------------------------------------------------------------------

function intv_width_TF = calcIntvWidthBelowThresh(idx0,idx1,intv_thresh,hilbEnergy,t)
% function returns the peak number that should be removed based on proximity and area under it
    intv_width_TF = [];

    intv_width = idx0(2:end) - idx1(1:end-1);
    
    % indices of peak intervals less than threshold specified
    pks_thresh = find(intv_width<intv_thresh);
    
    % if no peak intervals are identified, return an empty array
    if isempty(pks_thresh)
        intv_width_TF = [];
        return
    end
    
    % calculate differences between peak number to identify peaks that are right next to each other
    
    % calculate areas of peaks for comparison
    tmp_area = calcPeakArea(idx0,idx1,hilbEnergy,t);
            
    % if second peak is greater, choose to get rid of first psak 
    if tmp_area(pks_thresh(1)) < tmp_area(pks_thresh(1)+1)
        intv_width_TF = pks_thresh(1);
    else
        intv_width_TF = pks_thresh(1)+1;
    end

end

%% ----------------------------------------------------------------------------------------











