% running hrv analysis
numFiles = 6;
testnumber = [31
35
36
37
39
40
41
42
43
44
46
47
48
49
50];

m11 = zeros(length(testnumber),6);
m12 = zeros(length(testnumber),6);
m21 = zeros(length(testnumber),6);
m22 = zeros(length(testnumber),6);

for i=1:length(testnumber)
    disp(['Test number: ',num2str(testnumber(i))]);
    [timeDomainHRV1, timeDomainHRV2,timeDomainHRV3, timeDomainHRV4] = pre_processAudio2(testnumber(i), numFiles)
    
    % store values
    m11(i,:) = timeDomainHRV1;
    m12(i,:) = timeDomainHRV2;
    m21(i,:) = timeDomainHRV3;
    m22(i,:) = timeDomainHRV4;
end