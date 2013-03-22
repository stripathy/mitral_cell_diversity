% merge data from cluster split


numSplits = 4;
cnt = 1;
for i = 1:numSplits
    loadFile = ['split' num2str(i) 'Out.mat'];
    if exist(loadFile,'file') == 0
        continue;
    end
    load(loadFile);
    if cnt == 1
        inds = (cnt-1)*length(reconStructOut) + 1:cnt*length(reconStructOut);
%     reconStructInTemp(inds) = currReconStructIn;
        reconStructOutTemp(inds) = reconStructOut;
    else
        reconStructOutTemp(end+1:end+length(reconStructOut)) = reconStructOut;
    end
    cnt = cnt + 1;
end

clear reconStructOut; 

for i = 1:length(randInds)
    reconStructOut(i) = reconStructOutTemp(randInds==i);
end

% reconStructOut = reconStructOutTemp(randInds);
% % unrandomize indexes
% % [blah, randInds] = sort(rand(length(pop),1));
% for i = 1:length(randInds)
%     structInInds = (i-1)*numSignals+1:i*numSignals;
%     structOutInds = (randInds(i)-1)*numSignals+1:randInds(i)*numSignals;
%     reconStructOut(structOutInds) = reconStructOutTemp(structInInds);
% %     reconStructIn(structOutInds) = reconStructInTemp(structInInds);
% end

% reconStructOut = reconStructOutTemp;
clear reconStructTemp currReconStructIn reconStructOutTemp currInds reconStructInTemp i

