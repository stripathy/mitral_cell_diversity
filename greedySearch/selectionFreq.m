%% figure out the relative neuron frequencies

popSizeCon = 10; % largest pop size to consider

bestPopCntRep = zeros(44, size(bestPop,2));
for i = 1:size(bestPop,2)
    bestPopCntRep(:,i) = histc(bestPop(1:popSizeCon,i), 1:44);
    bestPopCntRepUni(:,i) = histc(unique(bestPop(1:popSizeCon,i)), 1:44);
end

bestPopFreqs = sum(bestPopCntRep,2);
bestPopFreqsUni = sum(bestPopCntRepUni,2);

figure;
bar(bestPopFreqs/numStimStats)
ylabel('average selection frequency');
xlabel('neuron number');

figure;
bar(bestPopFreqsUni/numStimStats)
ylabel('average selection frequency');
xlabel('neuron number');

bestPopFreqs

%% 

figure; hold on;
for i = 1:size(pMatrix,1)
    plot(pMatrix(i,1), pMatrix(i,2), '.', 'Color', colorCode(i,:,:), 'MarkerSize', 6 + 4*(bestPopFreqsUni(i)));
end

for i = 1:numSearchIters
    currNeuron = addedNeuron(i);
%     text(pMatrix(currNeuron,1)+.5, pMatrix(currNeuron,2), num2str(i));
end
xlabel('PCA 1');
ylabel('PCA 2');

%% compute neuron reliability - need to simulate first


%compute mean zero lag correlation between neuron spikes
corrBoxSize = 5;
numSignals = 400;
numSimRepeats = 10;
zeroLagCorrMatFull = zeros(numModels, numModels, numSignals);
parfor k = 1:numSignals
    %get spike trains for each neuron each trial
    spikeMat = [];
    for j = 1:length(cellModels)
        cellStimInd = find(simSpikesInds.stimInd == k & simSpikesInds.modelInd == j);
        currSpikes = simSpikes(cellStimInd).dat;
        spikeMat = vectCat(spikeMat, currSpikes);
    end
    tempCorrMat = myCorr1(spikeMat, corrBoxSize);
    [zeroLagCorrMatFull(:,:,k)] = corrMatProcess(tempCorrMat, numModels, numSimRepeats);
end
zeroLagCorrMat = nanmean(zeroLagCorrMatFull,3);
corrMat = zeroLagCorrMat;
symCorrMat = corrMat + triu(corrMat,1)';

for i = 1:numSignals/numSignalsPerStat
    inds = (1 + (i-1)*numSignalsPerStat):((i)*numSignalsPerStat);
    zeroLagCorrMat(:,:,i) = nanmean(zeroLagCorrMatFull(:,:,inds),3);
    cellReliab(:,i) = diag(zeroLagCorrMat(:,:,i));
end

