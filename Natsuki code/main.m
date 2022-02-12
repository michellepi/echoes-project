%% Main script
clear

%% Healthy challenge data classification results
% entire healthy folder path
Folder_He = '/Users/natsukib/OneDrive/3rd Year/Dissertation/Sound analysis/Heart sounds/Challenge healthy/';
[final_h] = final(Folder_He,1);

%% AS challenge data classification results
% entire as folder path
Folder_AS = '/Users/natsukib/OneDrive/3rd Year/Dissertation/Sound analysis/Heart sounds/Challenge aortic disease HS/';
[final_as] = final(Folder_AS,2);

%% SDR Classification threshold optimisation
% class matrix with sensitivity and specificity
class(:,1) = final_h';  % Spec
class(:,2) = final_as'; % Sens

% inverse specificity
for f = 1:length(class)
    inv_sp = 1 - class(:,1);
end

% plot ROC curve
figure
plot(inv_sp,class(:,2))
title('ROC Curve'); xlabel('1-Specificity'); ylabel('Sensitivity');
hold on

% plot reference line
hline = refline([1 0]); hline.Color = 'k'; hline.LineStyle = '-'; hline.HandleVisibility = 'off'; hold off;

% area under ROC curve
auc = abs(trapz(inv_sp(:),class(:,2)));
        