%% compare distn of best populations to an average distribution

% find the best pops of each stim distn
fileLists = dir('greedySearch*');
[fileLists] = sortGreedySearchNames(fileLists);

numStimStats = length(fileLists);
popCnt = 1;
popListBest = zeros(5, numStimStats*10);
for s = 1:numStimStats
    load(fileLists(s).name);
%     load(['greedySearch1',num2str(s)]);
    tempBestPop = getFinalPop(bestPopListAll);
    [bestPop(:,s)] = tempBestPop;
    
    for i = 1:length(tempBestPop)
        popListBest(1:i, popCnt) = sort(tempBestPop(1:i));
        popCnt = popCnt + 1;
    end
end

%% figure out how diverse each of the bestPops are
numStimStats = length(fileLists);

numBestPops = size(popListBest,2);
for i = 1:numBestPops
    tempPopList = nonzeros(popListBest(:,i));
    
    [stimDiff(:,i) histDiff(:,i) biasDiff(:,i)] = glmDiff(cellModels(tempPopList));
    [stimDiffUni(:,i) histDiffUni(:,i) biasDiffUni(:,i)] = glmDiff(cellModels(unique(tempPopList)));
end

bestAvgDivStim = mean(reshape(stimDiff,numSearchIters,numStimStats),2);
bestAvgDivHist = mean(reshape(histDiff,numSearchIters,numStimStats),2);
bestAvgDivBias = mean(reshape(biasDiff,numSearchIters,numStimStats),2);

stdErrDivStim = myStdErr(reshape(stimDiff,numSearchIters,numStimStats)');
stdErrDivHist = myStdErr(reshape(histDiff,numSearchIters,numStimStats)');
stdErrDivBias = myStdErr(reshape(biasDiff,numSearchIters,numStimStats)');

bestAvgDivStimUni = nanmean(reshape(stimDiffUni,numSearchIters,numStimStats),2);
bestAvgDivHistUni = nanmean(reshape(histDiffUni,numSearchIters,numStimStats),2);
bestAvgDivBiasUni = nanmean(reshape(biasDiffUni,numSearchIters,numStimStats),2);

stdErrDivStimUni = myStdErr(reshape(stimDiffUni,numSearchIters,numStimStats)');
stdErrDivHistUni = myStdErr(reshape(histDiffUni,numSearchIters,numStimStats)');
stdErrDivBiasUni = myStdErr(reshape(biasDiffUni,numSearchIters,numStimStats)');

%% figure out how diverse a randomly sampled pop is
numRandomPops = 500000;

numBestPops = size(popListBest,2);
parfor i = 1:numRandomPops
    popSize(i) = randi(numSearchIters);
    tempPopList = randsample(neuronsToSim, popSize(i), 1);    
    
    [stimDiffRand(:,i) histDiffRand(:,i) biasDiffRand(:,i)] = glmDiff(cellModels(tempPopList));
end

pctCI = 50;
pctCI2 = 90;
for i = 1:numSearchIters
%     popInds = find(
    meanDivStim(i) = mean(stimDiffRand(:,popSize == i));
    meanDivHist(i) = mean(histDiffRand(:,popSize == i));
    meanDivBias(i) = mean(biasDiffRand(:,popSize == i));
    
    meanStdDivStim(i) = std(stimDiffRand(:,popSize == i));
    meanStdDivHist(i) = std(histDiffRand(:,popSize == i));
    meanStdDivBias(i) = std(biasDiffRand(:,popSize == i));
      
    meanCIStim(i,:) = getConfInterval(stimDiffRand(:,popSize == i), pctCI);
    meanCIHist(i,:) = getConfInterval(histDiffRand(:,popSize == i), pctCI);
    meanCIBias(i,:) = getConfInterval(biasDiffRand(:,popSize == i), pctCI);
    
    meanCIStim2(i,:) = getConfInterval(stimDiffRand(:,popSize == i), pctCI2);
    meanCIHist2(i,:) = getConfInterval(histDiffRand(:,popSize == i), pctCI2);
    meanCIBias2(i,:) = getConfInterval(biasDiffRand(:,popSize == i), pctCI2);
end

%% find where observed mean best values are significantly less than random
% 
signifValsStim = find(bestAvgDivStim < meanCIStim2(:,1));
signifValsHist = find(bestAvgDivHist < meanCIHist2(:,1));
signifValsBias = find(bestAvgDivBias < meanCIBias2(:,1));

%% plot best diversity values vs randomly sampled ones

figure;
subplot(1,3,1); hold on;
% plot(1:numSearchIters, bestAvgDivStim, 1:numSearchIters, meanDivStim, 'g'); 
errorbar(1:numSearchIters, meanDivStim, meanCIStim(:,1)' - meanDivStim, meanCIStim(:,2)' - meanDivStim, 'g');
% errorbar(1:numSearchIters, meanDivStim, meanStdDivStim, 'g');
errorbar(1:numSearchIters, bestAvgDivStim, stdErrDivStim);
plot(2:numSearchIters, zeros(1, numSearchIters-1), 'r'); % represents homo diversity
plot(signifValsStim, 1.3, 'b*');
axis([1 numSearchIters+1 -.2 1.35]);
ylabel('Population diversity');
set(gca, 'XTick', [2 10])

subplot(1,3,2); hold on;
% plot(1:numSearchIters, bestAvgDivHist, 1:numSearchIters, meanDivHist, 'g'); 
errorbar(1:numSearchIters, meanDivHist, meanCIHist(:,1)' - meanDivHist, meanCIHist(:,2)' - meanDivHist, 'g');
% errorbar(1:numSearchIters, meanDivHist, meanStdDivHist, 'g');
errorbar(1:numSearchIters, bestAvgDivHist, stdErrDivHist);
plot(2:numSearchIters, zeros(1, numSearchIters-1), 'r'); % represents homo diversity
plot(signifValsHist, 45, 'b*');
axis([1 numSearchIters+1 -5 50]);
xlabel('Neuron addition number');
set(gca, 'XTick', [2 10])

subplot(1,3,3); hold on;
% plot(1:numSearchIters, bestAvgDivBias, 1:numSearchIters, meanDivBias, 'g'); 
errorbar(1:numSearchIters, meanDivBias, meanCIBias(:,1)' - meanDivBias, meanCIBias(:,2)' - meanDivBias, 'g');
% errorbar(1:numSearchIters, meanDivBias, meanStdDivBias, 'g');
errorbar(1:numSearchIters, bestAvgDivBias, stdErrDivBias);
% plot(signifValsBias, 18, 'b*');
plot(2:numSearchIters, zeros(1, numSearchIters-1), 'r'); % represents homo diversity
axis([1 numSearchIters+1 -3 25]);
set(gca, 'XTick', [2 10])


figure;
subplot(1,3,1);
errorbar(1:numSearchIters, meanDivStim, meanStdDivStim, 'g'); hold on;
% plot(1:numSearchIters, bestAvgDivStim, 1:numSearchIters, meanDivStim, 'g'); hold on;
errorbar(1:numSearchIters, bestAvgDivStim, stdErrDivStim);
errorbar(1:numSearchIters, bestAvgDivStimUni, stdErrDivStimUni, 'm');
axis([0 numSearchIters+1 0 1.25]);
ylabel('Population diversity');
set(gca, 'XTick', [1 10])

subplot(1,3,2);
errorbar(1:numSearchIters, meanDivHist, meanStdDivHist, 'g'); hold on;
% plot(1:numSearchIters, bestAvgDivHist, 1:numSearchIters, meanDivHist, 'g'); hold on;
errorbar(1:numSearchIters, bestAvgDivHist, stdErrDivHist);
errorbar(1:numSearchIters, bestAvgDivHistUni, stdErrDivHistUni, 'm');
axis([0 numSearchIters+1 0 40]);
xlabel('Neuron addition number');
set(gca, 'XTick', [1 10])

subplot(1,3,3);
errorbar(1:numSearchIters, meanDivBias, meanStdDivBias, 'g'); hold on;
% plot(1:numSearchIters, bestAvgDivBias, 1:numSearchIters, meanDivBias, 'g'); hold on;
errorbar(1:numSearchIters, bestAvgDivBias, stdErrDivBias);
errorbar(1:numSearchIters, bestAvgDivBiasUni, stdErrDivBiasUni, 'm');
axis([0 numSearchIters+1 0 20]);
set(gca, 'XTick', [1 10])


% compute a histogram over number of unique neurons in bestPops per stim



for i = 1:numStimStats
    numUniqueCells(i) = numel(unique(bestPop(:,i)));
end

histEdges = 0:10;
% figure;
[uniCnts, uniCenters] = histc(numUniqueCells,histEdges);
% bar(histEdges, uniCnts);

% count how many unique cells per pop one would expect by random chance
randPops = randi(length(neuronsToSim), 10, 20000);
parfor i = 1:20000
    numUniqueCellsRand(i) = numel(unique(randPops(:,i)));
end
[randCnts, randCenters] = histc(numUniqueCellsRand, 0:10);
uniqueCellsExpDistn = randCnts./sum(randCnts) * numStimStats;

hold on;
plot(0:10, uniqueCellsExpDistn, 'g')
xlabel('# unique neurons')
ylabel('count');

plotyy(histEdges, uniCnts, histEdges, uniqueCellsExpDistn, @bar, @plot)

% try to get error bars on expected number via random sampling
numRandSamples = 1000;
for i = 1:numRandSamples
    sampleInts = randsample(1:20000, 20000, 1);
    [randCnts, randCenters] = histc(numUniqueCellsRand(sampleInts), 0:10);
    uniqueCellsExpDistnResamp(i,:) = randCnts./sum(randCnts) * numStimStats;
end

for i = 1:1000
    for j = 1:numStimStats
        tempBestPop = randsample(44, 10, 1);
        randNumUnique(j) = numel(unique(tempBestPop));
    end
    tempHist(:,i) = histc(randNumUnique, histEdges);
end

expHist = mean(tempHist');
pctCI = 95;

for i = 1:length(histEdges)
    expCI(:,i) = getConfInterval(tempHist(i,:), 95);
end

signifVals = histEdges(uniCnts < expCI(1,:) | uniCnts > expCI(2,:));
% plot(histEdges, expHist)

figure; hold on;
bar(histEdges, uniCnts)
errorbar(histEdges, expHist, expCI(1,:) - expHist, expCI(2,:) - expHist, 'g');
plot(signifVals, 8, 'b*');
axis([0 11 0 9]);
xlabel('# unique neurons')
ylabel('count');
    




%% compute probability distn for monte carlo samples of diversity data

histEdges = linspace(min(stimDiffRand), max(stimDiffRand), 100);

probDistn = histc(stimDiffRand, histEdges)/numel(stimDiffRand);

meanVal = mean(stimDiffRand);
stdVal = std(stimDiffRand);
Y = normpdf(histEdges,meanVal,stdVal);
Y = Y/sum(Y);

probVals = normpdf(stimDiff(10:10:90), meanVal, stdVal)/sum(Y);

probValsRand = normpdf(randn(9,1)*stdVal + meanVal, meanVal, stdVal)/sum(Y);
sum(log(probValsRand))

mean(stimDiff(10:10:90))

for i = 1:10
%     probVals = normcdf(bestAvgDivStim', meanDivStim, meanStdDivStim);
    
    
    confIntStim(i,:) = norminv([.1], meanDivStim(i), meanStdDivStim(i));
end
