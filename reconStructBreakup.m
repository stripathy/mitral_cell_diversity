%% script to break up reconStructIn


%randomize indexes
[blah, randInds] = sort(rand(length(pop),1));
for i = 1:length(pop)
    structInInds = (i-1)*numSignals+1:i*numSignals;
    structOutInds = (randInds(i)-1)*numSignals+1:randInds(i)*numSignals;
    reconStructInSorted(structInInds) = reconStructIn(structOutInds);
end

numBreaks = 2;
numReconsPerBreak = length(reconStructInSorted)/numBreaks;

for i = 1:numBreaks
    currInds = (i-1)*numReconsPerBreak+1:i*numReconsPerBreak;
    currReconStructIn = reconStructInSorted(currInds);
    
    saveName = ['split' num2str(i)];
    save(saveName, 'cellModels', 'currReconStructIn', 'stimCovMat', 'outSignal');
    
end

%% break up for normApprox

%randomize indexes
[~, randInds] = sort(rand(length(reconStructIn),1));
reconStructInSorted = reconStructIn(randInds);

numBreaks = 4;
numReconsPerBreak = ceil(length(reconStructInSorted)/numBreaks);

for i = 1:numBreaks
    currInds = (i-1)*numReconsPerBreak+1:min(i*numReconsPerBreak, length(reconStructInSorted));
    currReconStructIn = reconStructInSorted(currInds);
    
    saveName = ['split' num2str(i)];
    save(saveName, 'cellModels', 'currReconStructIn', 'stimCovMat', 'outSignal');
    
end