function [sd_ratio,avg_sys_amp,avg_dias_amp,count] = AmplitudeMetrics(ampt_12,ampt_21,ratio)

%% Amplitudes
ampt_12 = ampt_12(ampt_12>0); ampt_21 = ampt_21(ampt_21>0);
avg_sys_amp = mean(ampt_12);
avg_dias_amp = mean(ampt_21);

%% making sure everything is the right size
diff1_ = length(ampt_12) - length(ampt_21);
if length(ampt_12) > length(ampt_21)
    ampt_21 = padarray(ampt_21,diff1_,0,'post');
elseif length(ampt_21) > length(ampt_12)
    ampt_12 = padarray(ampt_12,-diff1_,0,'post');
end

%% number of systoles that have higher amplitudes than their corresponding diastoles
count = 0;
for p = 1:length(ampt_12)
    if ampt_12(p) > 0 || ampt_21(p) > 0
        if ampt_12(p) > ampt_21(p) && (ampt_12(p) - ampt_21(p)) > 0.75 * ampt_21(p)
            count = count + 1;
        end
    end
end

%% Amplitude metric used
% every heart cycle's systolic/diastolic ratio in each recording
ratioTot = ratio;
ratioTot(ratioTot==0) = NaN;
sd_ratio = mean(ratioTot,'omitnan');

end

