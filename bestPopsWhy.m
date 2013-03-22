% why are some populations the best?

%% what differentiates the good from the other populations?


% bestPopsThresh = 10; % the threshold to look for "good pops"

%find the bestest pops
lookList = [1 3];
bestPopsAllStims = find(histc(nonzeros(bestPops(lookList,1:bestPopsThresh)),0:1:numPops+1) >= numel(lookList))-1;

%sort each bestest pop in order of lowest sum position across each stim
tempSum= zeros(size(bestPopsAllStims,1), 1);
for j = 1:size(bestPopsAllStims,1)
    for i = 1:numStimStats
        tempSum(j) = tempSum(j) + find(bestPops(i,:) == bestPopsAllStims(j));
    end
end
[sortedTempSum sortedTempSumInds] = sort(tempSum);
bestPopList = popList(:,bestPopsAllStims(sortedTempSumInds));

% count which neurons occur in bestPopList
histEdges = 1:numModels;
[bestPopsAllStimsHist] = histc(nonzeros(bestPopList(:,:)),histEdges);
% figure; 
% bar(1:length(neuronsToSim), count(1:end));
% xlabel('neuron index'); ylabel('counts');

%find the bestest pops for each stim
for i = 1:numStimStats
    bestNeuronHist(i,:) = histc(nonzeros(popList(:,bestPops(i,1:bestPopsThresh))), histEdges);
    bestNeuronHist(i,:) = bestNeuronHist(i,:)/sum(bestNeuronHist(i,:));
end

% plot a histogram over how often each neuron is in the best population for
% each stim stat

bestNeuronFreq = bestPopsAllStimsHist/sum(bestPopsAllStimsHist);
figure;
subplot(221)
% plot(histEdges, bestNeuronHist)
hold on;
plot(histEdges, bestNeuronFreq, 'k' ,'LineWidth', 2);
xlabel('neuron identity'); ylabel('probability');
axis([0 45 0 max(bestNeuronFreq)+.02]);

figure;
hold on;
plot(histEdges, bestNeuronHist(1, :), 'r')
plot(histEdges, bestNeuronHist(3, :), 'b')
xlabel('neuron identity'); ylabel('probability');
legend('high freq', 'low freq');

% figure out 2nd ord relationships among best neurons - just count
pairwiseFreqMat = zeros(numModels, numModels);
for i = 1:size(bestPopList,2)
    tempCombs = combnk(bestPopList(:,i), 2);
%     tempCombs = combnk(1:5, 2);
    for j = 1:size(tempCombs,1)
        pairwiseFreqMat(tempCombs(j,1), tempCombs(j,2)) = pairwiseFreqMat(tempCombs(j,1), tempCombs(j,2)) + 1;
    end
end
pairwiseFreqMat = pairwiseFreqMat/sum(nonzeros(pairwiseFreqMat));

%% is the pariwise histo diff from the histo via 1st order assumptions?

pairwiseFreqMat1st = bestNeuronFreq * bestNeuronFreq';
pairwiseFreqMat1st = tril(pairwiseFreqMat1st, -1)' + pairwiseFreqMat1st;
pairwiseFreqMat1st = pairwiseFreqMat1st - tril(pairwiseFreqMat1st, -1);

diffFreqMat = pairwiseFreqMat - pairwiseFreqMat1st;
figure; imagesc(diffFreqMat); colorbar;

rand1stSamples = randsample(numModels, 50000, 'true', bestNeuronFreq);

rand1stSamples = sort(reshape(rand1stSamples, 5, 10000));

% figure out 2nd ord relationships among best neurons - just count
pairwiseFreqMatSamp = zeros(numModels, numModels);
for i = 1:size(rand1stSamples,2)
    tempCombs = combnk(rand1stSamples(:,i), 2);
%     tempCombs = combnk(1:5, 2);
    for j = 1:size(tempCombs,1)
        pairwiseFreqMatSamp(tempCombs(j,1), tempCombs(j,2)) = pairwiseFreqMatSamp(tempCombs(j,1), tempCombs(j,2)) + 1;
    end
end
pairwiseFreqMatSamp = pairwiseFreqMatSamp/sum(nonzeros(pairwiseFreqMatSamp));

%% how diverse is each of the good pops relative to the full distn?

for i = 1:numPops
    popDiversity(i) = numel(unique(popList(:,i)))/maxPopSize;
end

% compute mean firing rates of each neuron
for j = 1:length(cellModels)
    fullSpikeMat = [];
    for k = 1:numSignals
        cellStimInd = k + (j-1)*numSignals;
        currSpikes = simSpikes(cellStimInd).dat;
%         [blah1, blah1, tempReliab(k)] = gammaCoincFact(currSpikes, currSpikes, testInterval, deltaLen);
        fullSpikeMat = vectCat(fullSpikeMat, currSpikes);
    end
    
    for i = 1:numStimStats
        stimInds = 1+(i-1)*numSignalsPerStat*numSimRepeats:(i)*numSignalsPerStat*numSimRepeats;
        meanFR(j,i) = nnz(fullSpikeMat(:,stimInds))/(numSimRepeats*numSignalsPerStat*slen*.001);
    end
    
%     cellReliab(j)  = mean(tempReliab);
end
meanFRPerStim = meanFR;
meanFR = mean(meanFR,2);

%compute mean zero lag correlation between neuron spikes
corrBoxSize = 5;
zeroLagCorrMatFull = zeros(numModels, numModels, numSignals);
parfor k = 1:numSignals
    %get spike trains for each neuron each trial
    spikeMat = [];
    for j = 1:length(cellModels)
        cellStimInd = k + (j-1)*numSignals;
        currSpikes = simSpikes(cellStimInd).dat;
        spikeMat = vectCat(spikeMat, currSpikes);
    end
    tempCorrMat = myCorr1(spikeMat, corrBoxSize);
    [zeroLagCorrMatFull(:,:,k)] = corrMatProcess(tempCorrMat, numModels, numSimRepeats);
end

for i = 1:numStimStats
    inds = 1+(i-1)*numSignalsPerStat:(i)*numSignalsPerStat;
    zeroLagCorrMat(:,:,i) = mean(zeroLagCorrMatFull(:,:,inds),3);
end
corrMat = zeroLagCorrMat;
meanCorrMat = mean(corrMat,3);

% use only hetero pops
% popList = popList(:,45:end);
% numPops = 894 - 44;
% bestPopsAllStims = bestPopsAllStims - 44;

% whats the mean FR for each population?
for i = 1:numPops
    popMeanFR(i) = mean(meanFR(popList(:,i)));
end

% what's the mean pairwise correlation for each pop?
for i = 1:numPops
    tempPerms = combnk(popList(:,i),2);
    for j = 1:length(tempPerms)
        popMeanCorrAll(j,i) = meanCorrMat(tempPerms(j,1), tempPerms(j,2));
    end
end
popMeanCorr = mean(popMeanCorrAll);

% mean reliability of each neuron in the pop?
meanReliab = diag(meanCorrMat);
for i = 1:numPops
    popMeanReliab(i) = mean(meanReliab(popList(:,i)));
end

% plot histograms of mean quantities vs those for best pops

%mean FR histo
histEdges = linspace(5, 60, 20);
[meanFRHist] = histc(popMeanFR,histEdges)/length(popMeanFR);
[bestPopsMeanFRHist] = histc(popMeanFR(bestPopsAllStims), histEdges)/length(popMeanFR(bestPopsAllStims));
% figure;
subplot(222)
plot(histEdges, [meanFRHist; bestPopsMeanFRHist]);
xlabel('Mean firing rate (Hz)'); ylabel('probability');
axis tight;

%mean reliab histo
histVals = popMeanReliab;
histEdges = linspace(min(histVals), max(histVals), 15);
[fullHist] = histc(histVals,histEdges)/length(histVals);
[bestPopsHist] = histc(histVals(bestPopsAllStims), histEdges)/length(histVals(bestPopsAllStims));

% figure;
subplot(223)
plot(histEdges, [fullHist; bestPopsHist]);
xlabel('Mean reliability'); ylabel('probability');
axis tight;

%mean pairwise corr histo
histEdges = linspace(min(popMeanCorr), max(popMeanCorr), 15);
[meanCorrHist] = histc(popMeanCorr,histEdges)/length(popMeanCorr);
[bestPopsMeanCorrHist] = histc(popMeanCorr(bestPopsAllStims), histEdges)/length(popMeanCorr(bestPopsAllStims));

subplot(224)
plot(histEdges, [meanCorrHist; bestPopsMeanCorrHist]);
xlabel('Mean pairise corr'); ylabel('probability');
axis tight;

% what's the mean difference among glm parameters of each pop?
stimDiffMat = zeros(numModels, numModels);
for i = 1:numModels
    for j = i+1:numModels
        [stimDiffMat(i,j) histDiffMat(i,j) biasDiffMat(i,j)] = glmDiff(cellModels(i), cellModels(j));
    end
end

for i = 1:numPops
    tempPerms = combnk(popList(:,i),2);
    for j = 1:length(tempPerms)
        [stimDiff(j,i) histDiff(j,i) biasDiff(j,i)] = glmDiff(cellModels(tempPerms(j,1)), cellModels(tempPerms(j,2)));
    end
end
meanStimDiff = mean(stimDiff); 
meanHistDiff = mean(histDiff); 
meanBiasDiff = mean(biasDiff); 
stdStimDiff = std(stimDiff);

%mean glm diffs histograms
histVals = meanStimDiff;
histEdges = linspace(min(histVals), max(histVals), 20);
[fullHist] = histc(histVals,histEdges)/length(histVals);
[bestPopsHist] = histc(histVals(bestPopsAllStims), histEdges)/length(histVals(bestPopsAllStims));

figure;
subplot(221)
plot(histEdges, [fullHist; bestPopsHist]);
xlabel('Pop mean stim diff'); ylabel('probability');


histVals = meanHistDiff;
histEdges = linspace(min(histVals), max(histVals), 20);
[fullHist] = histc(histVals,histEdges)/length(histVals);
[bestPopsHist] = histc(histVals(bestPopsAllStims), histEdges)/length(histVals(bestPopsAllStims));

subplot(222)
plot(histEdges, [fullHist; bestPopsHist]);
xlabel('Pop mean hist diff'); ylabel('probability');

histVals = meanBiasDiff;
histEdges = linspace(min(histVals), max(histVals), 20);
[fullHist] = histc(histVals,histEdges)/length(histVals);
[bestPopsHist] = histc(histVals(bestPopsAllStims), histEdges)/length(histVals(bestPopsAllStims));

subplot(223)
plot(histEdges, [fullHist; bestPopsHist]);
xlabel('Pop mean bias diff'); ylabel('probability');

%mean population single cell infos
for i = 1:numPops
    popInfo(:,i) = mean(totalEntPerSec(:,popList(:,i)),2);
end
popInfoMean = mean(popInfo);


histVals = popInfoMean;
histEdges = linspace(min(histVals), max(histVals), 20);
[fullHist] = histc(histVals,histEdges);%/length(histVals);
[bestPopsHist] = histc(histVals(bestPopsAllStims), histEdges);%/length(histVals(bestPopsAllStims));

figure;
plot(histEdges, [fullHist; bestPopsHist]);
xlabel('Pop mean single cell info'); ylabel('probability');

% plot GLM parms for example pops
tempPopInd = 822;
currCellInds = popList(:,tempPopInd);
% currCellInds = [10 27 33 41 44];
% currCellInds = [27 33];
colorCode = colormap(jet(length(keepCells)));

figure; 
allCellInds = unique(nonzeros(currCellInds));
for i = 1:length(allCellInds)
    
    subplot(1,5,1:2);plot(-49:0, cellModels(allCellInds(i)).k, 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 2); axis([-50 0 -.5 1.2]); hold on;
    subplot(1,5,3:4);plot(cellModels(1).iht,log10(exp(cellModels(allCellInds(i)).ihbas*cellModels(allCellInds(i)).ih)), 'Color', colorCode(allCellInds(i),:), 'LineWidth' , 2); hold on;
    axis ([0 60 -1.6 2.2]);
    yTickPos = [-1:2]; yTickLabels = [.1 1 10 100];
    set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabels);
    subplot(1,5,5);plot(1, log10(exp(cellModels(allCellInds(i)).dc)),  '*', 'Color', colorCode(allCellInds(i),:)); hold on;
    axis ([0 2 -1.6 2.2]);
    set(gca, 'YTick', yTickPos, 'YTickLabel', [], 'XTick', []);
end

%% try regressing out effects

% see how stim 1 depends on predictor variables
for i = 1:numPops
    popMeanFR(i) = mean(meanFRPerStim(popList(:,i),1));
end

plot(popMeanFR, reconRMSEMean(1,:),'.')

X = [ones(numPops, 1) popMeanFR'];
y = reconRMSEMean(1,:)';
betahat = inv(X'*X)*X'*y;
resids = y - X*betahat;
plot(popMeanFR, reconRMSEMean(1,:),'.' )

pctFR = prctile(popMeanFR, 90);
hold on;
plot([1 1]*pctFR, [.5 .9], 'g');

bestPopInds = bestPops(1,1:100);
plot(popMeanFR(bestPopInds), reconRMSEMean(1,bestPopInds),'r.' )

stimInd = 1;
for i = 1:numPops
    popMeanFR(i) = mean(meanFRPerStim(popList(:,i),stimInd));
end

indepVar = popMeanReliab;
depVar = reconRMSEMean(stimInd,:);
% depVar = resids;
pctVar = prctile(indepVar, 90);
bestPopInds = bestPops(stimInd,1:100);
bestPopIndsBoth = intersect(bestPops(1,1:100), bestPops(3,1:100));
pctCorrect = sum(indepVar(bestPopInds)>=pctVar)/numel(bestPopInds);

figure;
plot(indepVar, depVar,'.' )
hold on;
%draw seperation line 
plot([1 1]*pctVar, [min(depVar) max(depVar)], 'g');
%draw top x % pops
plot(indepVar(bestPopInds), depVar(bestPopInds),'r.' )

% try plotting two indep vars and class category by color
indepVar1 = popMeanFR;
indepVar2 = meanStimDiff;

figure;
plot(indepVar1, indepVar2, 'b.');
hold on;
plot(indepVar1(bestPopInds), indepVar2(bestPopInds), 'r.');

% % color cells directly by FR
colorCodeTotal = colormap(hot(max(floor(depVar*1000)) - min(floor(depVar*1000)) + 1));
for i = 1:numPops
    colorCode(i,:) = colorCodeTotal(floor(depVar(i)*1000) - min(floor(depVar*1000)) +1,:);
end

figure;
hold on;
for i = 1:numPops
    plot(indepVar1(i), indepVar2(i), '.', 'Color', colorCode(i,:));
end
colormap hot; caxis([min(depVar) max(depVar)]);

color
% plot(indepVar1(bestPopIndsBoth), indepVar2(bestPopIndsBoth), 'g.');

plot(popMeanFR, resids,'.');


plot(popMeanReliab, reconRMSEMean(1,:),'.')

% try regressing and looking at r^2 values

stimInd = 1;
for i = 1:numPops
    popMeanFR(i) = mean(meanFRPerStim(popList(:,i),stimInd));
    tempPerms = combnk(popList(:,i),2);
    for j = 1:length(tempPerms)
        popMeanCorrAll(j,i) = corrMat(tempPerms(j,1), tempPerms(j,2), stimInd);
    end    
end
popMeanCorr = mean(popMeanCorrAll);

meanReliab = diag(corrMat(:,:, stimInd));
for i = 1:numPops
    popMeanReliab(i) = mean(meanReliab(popList(:,i)));
end

popMakeup = histc(popList,1:44);

% indepVar = popMeanReliab;
depVar = reconRMSEMean(stimInd,:);
% depVar = mean(reconRMSEMean([1 3],:));
indepVar1 = popMeanFR;
indepVar2 = popMeanReliab;
% indepVar3 = popInfo(2,:);
indepVar4 = popMeanCorr;
indepVar5 = meanHistDiff;

y = depVar';
X = [ones(1, 894)' popMeanFR' zscore(meanHistDiff')];

[b, bint, r, rint,  stats]= regress(y, X);
stats(1)

plot(popMeanReliab, reconRMSEMean(stimInd,:),'.')

stim1PctVarExp = [69 72 86 85 77 79 94];
stim3PctVarExp = [50 60 66 62 63 68 89];
labels = {'FR', 'FR+meanCorr', 'FR+stimDiff+histDiff', 'FR+stimDiff',...
    'FR+histDiff', 'popMakeup', 'popMakeup+meanCorr'};
figure;bar(stim3PctVarExp);
set(gca,'XTickLabel', labels);

%% screw PCA space, just avg filters and compare
clear popMeanParms* popVarParms*
for i = 1:numPops
    popMeanParmsK(:,i) = mean([kMat(:,popList(:,i))']',2);
    popVarParmsK(:,i) = var([kMat(:,popList(:,i))']', [], 2);
    popMeanParmsH(:,i) = mean([ihMat(:,popList(:,i))']',2);
    popVarParmsH(:,i) = var([ihMat(:,popList(:,i))']', [], 2);
    popMeanParmsB(:,i) = mean([dcComp(:,popList(:,i))']',2);
    popVarParmsB(:,i) = var([dcComp(:,popList(:,i))']', [], 2);
%     popMeanParms(:,i) = mean([kMat(:,popList(:,i))' ihMat(:,popList(:,i))' dcComp(:,popList(:,i))']',2);
%     popVarParms(:,i) = var(projFilters(popList(:,i),:));
%     popCovParms(:,:,i) = cov(projFilters(popList(:,i),:));
end

% plot(1:50, mean(popMeanParms(:, bestPops(3,1:10)),2))

numBest = 80;
figure;
subplot(1,5,1:2);
hold on;
meanParms = mean(popMeanParmsK(:, bestPops(1,1:numBest)),2);
stdParms = mean(sqrt(popVarParmsK(:, bestPops(1,1:numBest))),2);
plot(tptsK, meanParms, 'r', 'LineWidth', 2);
plot(tptsK, repmat(meanParms,1,2)' + [-1 1]'*stdParms','r--')

meanParms = mean(popMeanParmsK(:, bestPops(2,1:numBest)),2);
stdParms = mean(sqrt(popVarParmsK(:, bestPops(2,1:numBest))),2);
plot(tptsK, meanParms, 'b', 'LineWidth', 2);
plot(tptsK, repmat(meanParms,1,2)' + [-1 1]'*stdParms','b--')
axis([-25 0 -.2, 1]);
set(gca, 'XTick', -25:5:0);
xlabel('Time (ms)');
ylabel('Stim units (a.u.)');

subplot(1,5,3:4);
hold on;
meanParms = mean(popMeanParmsH(:, bestPops(1,1:numBest)),2);
stdParms = mean(sqrt(popVarParmsH(:, bestPops(1,1:numBest))),2);
plot(ggAll(1).iht, meanParms, 'r', 'LineWidth', 2);
plot(ggAll(1).iht, repmat(meanParms,1,2)' + [-1 1]'*stdParms','r--')

meanParms = mean(popMeanParmsH(:, bestPops(2,1:numBest)),2);
stdParms = mean(sqrt(popVarParmsH(:, bestPops(2,1:numBest))),2);
plot(ggAll(1).iht, meanParms, 'b', 'LineWidth', 2);
plot(ggAll(1).iht, repmat(meanParms,1,2)' + [-1 1]'*stdParms','b--')
axis ([0 60 -1.6 2.2]);
yTickPos = [-1:2]; yTickLabels = [.1 1 10 100];
set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabels);
set(gca, 'XTick', 0:10:60);
xlabel('Time (ms)');
ylabel('Gain');

subplot(1,5,5);
hold on;
meanParms = mean(popMeanParmsB(:, bestPops(2,1:numBest)),2);
stdParms = mean(sqrt(popVarParmsB(:, bestPops(2,1:numBest))),2);
errorbar(1, meanParms, stdParms, 'xb', 'MarkerSize', 8);
meanParms = mean(popMeanParmsB(:, bestPops(1,1:numBest)),2);
stdParms = mean(sqrt(popVarParmsB(:, bestPops(1,1:numBest))),2);
errorbar(1, meanParms, stdParms, 'xr', 'MarkerSize', 8);
% plot(1, meanParms, 'xr',  1, repmat(meanParms,1,2)' + [-1 1]'*stdParms','r--')
axis ([.8 1.2 -1.6 2.2]);
set(gca, 'YTick', yTickPos, 'YTickLabel', []);
set(gca, 'XTick', []);


% meanParms = mean(popMeanParms(:, bestPops(4,1:numBest)),2);
% stdParms = mean(sqrt(popVarParms(:, bestPops(4,1:numBest))),2);
% plot(1:length(meanParms), meanParms, 'xg', 1:length(meanParms), repmat(meanParms,1,2)' + [-1 1]'*stdParms','g--')



figure;
imagesc([ mean(popMeanParms(:,bestPops(1,1:10)),2) mean(popVarParms(:,bestPops(1,1:10)),2)]);
colorbar; caxis([-.3 1.5]);

figure;
imagesc([ mean(popMeanParms(:,bestPops(3,1:10)),2) mean(popVarParms(:,bestPops(3,1:10)),2)]);
colorbar; caxis([-.3 1.5]);

[vals bestPopIndsBoth] = sort(mean(reconRMSEMean([1 3],:)));

figure;
imagesc([ mean(popMeanParms(:,bestPopIndsBoth(1:10)),2) mean(popVarParms(:,bestPopIndsBoth(1:10)),2)]);
colorbar; caxis([-.3 1.5]);

%plot mean glm parameter values
bestParmsMean = mean(popMeanParms(:,bestPops(1,1:10)),2);
bestParmsVar = mean(popVarParms(:,bestPops(1,1:10)),2);

% bestParmsMean = mean(popMeanParms(:,bestPopIndsBoth(1:10)),2);
% bestParmsVar = mean(popVarParms(:,bestPopIndsBoth(1:10)),2);


meanParmVals = bestParmsMean.*(projFilterStds)';
meanParmFilts = meanParmVals(1:numEigsK)'*outputPCsK(:,1:numEigsK)';

stdParmVals = sqrt(bestParmsVar).*(projFilterStds)';
stdParmFilts = repmat(meanParmFilts, 2, 1)  + [-1 1]'*stdParmVals(1:numEigsK)'*outputPCsK(:,1:numEigsK)';

figure;
plot(1:50, meanParmFilts, 1:50, stdParmFilts, '--');


meanParmVals = bestParmsMean.*(projFilterStds)';
meanParmFilts = meanParmVals(numEigsK+1:numEigsK + numEigsH)'*outputPCsH(:,1:numEigsH)';

stdParmVals = sqrt(bestParmsVar).*(projFilterStds)';
stdParmFilts = repmat(meanParmFilts, 2, 1)  + [-1 1]'*stdParmVals(numEigsK+1:numEigsK + numEigsH)'*outputPCsH(:,1:numEigsH)';

figure;
plot(1:531, meanParmFilts, 1:531, stdParmFilts, '--');


figure;
plot(1:50, mean(mean(nonzeros(popList(:,bestPops(1,1:10))))));
for i = 1:10
    avgStimFilt = avgStimFilt + popList

% plotting of means and covariance matrices
figure;
imagesc([ mean(popMeanParms(:,bestPops(3,1:10)),2) zeros(6,1) mean(popCovParms(:,:,bestPops(3,1:10)),3)]);
colorbar;

figure;
imagesc([ mean(popMeanParms(:,bestPops(1,1:10)),2) zeros(6,1) mean(popCovParms(:,:,bestPops(1,1:10)),3)]);
colorbar;

figure;
imagesc([ mean(popMeanParms(:,bestPops(3,1:10)),2) zeros(6,1) mean(popCovParms(:,:,bestPops(3,1:10)),3)]);
colorbar;

%% in reduced GLM parm space, is there an optimal set of filters per stim?

for i = 1:numPops
    popMeanParms(:,i) = mean(projFilters(popList(:,i),:));
    popVarParms(:,i) = var(projFilters(popList(:,i),:));
    popCovParms(:,:,i) = cov(projFilters(popList(:,i),:));
end

figure;
imagesc([ mean(popMeanParms(:,bestPops(1,1:10)),2) mean(popVarParms(:,bestPops(1,1:10)),2)]);
colorbar; caxis([-.3 1.5]);

figure;
imagesc([ mean(popMeanParms(:,bestPops(3,1:10)),2) mean(popVarParms(:,bestPops(3,1:10)),2)]);
colorbar; caxis([-.3 1.5]);

[vals bestPopIndsBoth] = sort(mean(reconRMSEMean([1 3],:)));

figure;
imagesc([ mean(popMeanParms(:,bestPopIndsBoth(1:10)),2) mean(popVarParms(:,bestPopIndsBoth(1:10)),2)]);
colorbar; caxis([-.3 1.5]);

%plot mean glm parameter values
bestParmsMean = mean(popMeanParms(:,bestPops(1,1:40)),2);
bestParmsVar = mean(popVarParms(:,bestPops(1,1:40)),2);

% bestParmsMean = mean(popMeanParms(:,bestPopIndsBoth(1:10)),2);
% bestParmsVar = mean(popVarParms(:,bestPopIndsBoth(1:10)),2);


meanParmVals = bestParmsMean.*(projFilterStds)';
meanParmFilts = meanParmVals(1:numEigsK)'*outputPCsK(:,1:numEigsK)';

stdParmVals = sqrt(bestParmsVar).*(projFilterStds)';
stdParmFilts = repmat(meanParmFilts, 2, 1)  + [-1 1]'*stdParmVals(1:numEigsK)'*outputPCsK(:,1:numEigsK)';

figure;
plot(1:50, meanParmFilts, 1:50, stdParmFilts, '--');


meanParmVals = bestParmsMean.*(projFilterStds)';
meanParmFilts = meanParmVals(numEigsK+1:numEigsK + numEigsH)'*outputPCsH(:,1:numEigsH)';

stdParmVals = sqrt(bestParmsVar).*(projFilterStds)';
stdParmFilts = repmat(meanParmFilts, 2, 1)  + [-1 1]'*stdParmVals(numEigsK+1:numEigsK + numEigsH)'*outputPCsH(:,1:numEigsH)';

figure;
plot(1:531, meanParmFilts, 1:531, stdParmFilts, '--');
%% how informative is each neuron on each stim?
% I actually need to compute this for the other stims, for now use infos
% from old stuff

figure;
plot(1:44, totalEntPerSec(1:44), '.')

figure;
plot(bestNeuronFreq, totalEntPerSec(1:44), '.')

colorCode = colormap(jet(length(keepCells)));
figure;
hold on;
for i = 1:numModels
    plot(totalEntPerSec(1,i), totalEntPerSec(2,i),  '*', 'Color', colorCode(i,:));
    text(totalEntPerSec(1,i)+1, totalEntPerSec(2,i)+1,num2str(i), 'FontSize', 6);
end
xlabel('Info per neuron: stim 1 (bits/sec)');
ylabel('Info per neuron: stim 3 (bits/sec)');


% mean info of the neurons in the pop?
meanPopInfo = diag(meanCorrMat);
for i = 1:numPops
%     meanPopInfo(i) = mean(mean(totalEntPerSec(:,popList(:,i))));
    meanPopInfo(i) = mean(totalEntPerSec(1,popList(:,i)));
end

%mean reliab histo
histVals = meanPopInfo;
histEdges = linspace(min(histVals), max(histVals), 15);
[fullHist] = histc(histVals,histEdges)/length(histVals);
[bestPopsHist] = histc(histVals(bestPopsAllStims), histEdges)/length(histVals(bestPopsAllStims));

figure;
% subplot(223)
plot(histEdges, [fullHist bestPopsHist]);
xlabel('Mean single cell info'); ylabel('probability');
axis tight;

[vals, sorted] = sort(totalEntPerSec(1,:),2);
top5PopHigh = sorted(end-4:end);
[vals, sorted] = sort(totalEntPerSec(2,:),2);
top5PopLow = sorted(end-4:end);

entStimStats = [1 3];
for i = 1:length(entStimStats)
    totalEntPerSecNorm(i,:) = totalEntPerSec(i,:)/max(totalEntPerSec(i,:),[],2);
end
[vals, sorted] = sort(mean(totalEntPerSecNorm),2);
top5PopBoth = sorted(end-4:end);

[vals, sorted] = sort(mean(totalEntPerSec),2);
top5PopBothNotNorm = sorted(end-4:end);

