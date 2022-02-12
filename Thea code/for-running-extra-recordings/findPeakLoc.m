function [idx0 idx1 hilbEnergy_segment] = findPeakLoc(hilbEnergy_segment)
% find locations of peaks
    ne0 = find(hilbEnergy_segment~=0);                               % Find non-zero elements
    idx0 = unique([ne0(1); ne0(diff([0; ne0])>1)]);                         % Start indices of non-zero elements
    end_idxs = find(diff([0; ne0])>1)-1;
    
    if ne0(1) ~= 1
        end_idxs(1) = [];
    end
    
    idx1 = ne0([end_idxs; length(ne0)]);                                    % End indices of non-zero elements

    
    % check if any peaks are single-valued
    for i = 1:length(idx0)
        if idx0(i) == idx1(i)
            % remove single-valued peak from indices and segment array
            hilbEnergy_segment(idx0(i)) = 0;
            idx0(i) = 0; idx1(i) = 0;
        end
    end
    idx0(idx0==0) = []; idx1(idx1==0) = [];
    
end