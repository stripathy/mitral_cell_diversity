% make figures for greedySearch stuff, for 5x more neurons

%% load the folder containing data
cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\greedySearch\finSims\sim2');

% load data from greedySearch iterations across stimuli
% these get stored in bestPop

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

%% for an example stim stat, plot the iterative ordering of the neurons

% first compute pca space of neurons
% load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\drugAllFitsFinalggAll.mat';

cellModels = ggAll;
keepCells = 1:length(ggAll);
numModels = length(ggAll);

clear kMat ihMat
histInds = 101:length(ggAll(1).ihbas);
for j = 1:length(keepCells)
    currCellIds(j) = keepCells(j);%find(allGroupInds(i,j)==cellId);
%     staMat(:,j) = ggAll(currCellIds(j)).sta;
    kMat(:,j) = ggAll(currCellIds(j)).k;
    temp = log10(exp(ggAll(currCellIds(j)).ihbas*ggAll(currCellIds(j)).ih));
    ihMat(:,j) = temp(histInds);
%     ihMat(:,j) = exp(ggAll(currCellIds(j)).ihbas*ggAll(currCellIds(j)).ih);
    dcComp(j) = log10(exp(ggAll(currCellIds(j)).dc));
end

allFilts = [kMat' downsample(ihMat,10)' dcComp']';
allFilts = zscore(allFilts);

data = allFilts;
allFiltsMean = mean(allFilts,2);
pcaEigs = 2;
[outputPCsA,pMatrixA, rankEigen] = princomp(data');
pctVarExpA = cumsum(rankEigen)./sum(rankEigen);

% plot neurons along first 2 PCs
exampleStimStats = [2 4];
neuronDistnPan = panel();
neuronDistnPan.pack('h', length(exampleStimStats));

for i = 1:length(exampleStimStats)

currStimStat = exampleStimStats(i);
currBestPop = bestPop(:, currStimStat);

pMatrix = pMatrixA;
outputPCs = outputPCsA;
colorCode = colormap(jet(length(keepCells)));
neuronDistnPan(i).select();

hold on;
for i = 1:size(pMatrix,1)
    plot(pMatrix(i,1) + randn*(.1), pMatrix(i,2)+ randn*(.1), '.', 'Color', colorCode(i,:,:), 'MarkerSize', 10);
end
% scatter(pMatrix(:,1), pMatrix(:,2), sizeVec, colorCode, 'filled');
for i = 1:numSearchIters
    currNeuron = currBestPop(i);
    text(pMatrix(currNeuron,1)+.5, pMatrix(currNeuron,2), num2str(i), 'FontSize', 6);
end
xlabel('PCA 1');
ylabel('PCA 2');
axis([-10 10 -4 6]);

end
neuronDistnPan.fontsize = 6;
%% figure out how diverse each of the bestPops are
removeStims = [7 9];
includeStims = setdiff(1:numStimStats, removeStims);
numStimStats = length(fileLists);

numBestPops = size(popListBest,2);
clear *Diff;
for i = 1:numBestPops
    tempPopList = nonzeros(popListBest(:,i));
    
    [stimDiff(:,i) histDiff(:,i) biasDiff(:,i)] = glmDiff(cellModels(tempPopList));
%     [stimDiffUni(:,i) histDiffUni(:,i) biasDiffUni(:,i)] = glmDiff(cellModels(unique(tempPopList)));
end

stimDiffAll = reshape(stimDiff,numSearchIters,numStimStats)';
histDiffAll = reshape(histDiff,numSearchIters,numStimStats)';
biasDiffAll = reshape(biasDiff,numSearchIters,numStimStats)';

stimDiffAll = stimDiffAll(includeStims,:);
histDiffAll = histDiffAll(includeStims,:);
biasDiffAll = biasDiffAll(includeStims,:);

bestAvgDivStim = mean(stimDiffAll,1);
bestAvgDivHist = mean(histDiffAll,1);
bestAvgDivBias = mean(biasDiffAll,1);

stdErrDivStim = myStdErr(stimDiffAll);
stdErrDivHist = myStdErr(histDiffAll);
stdErrDivBias = myStdErr(biasDiffAll);


for i = 1:numSearchIters
    
bestCIStim(i,:) = getConfInterval(stimDiffAll(:,i), pctCI);
bestCIHist(i,:) = getConfInterval(histDiffAll(:,i), pctCI);
bestCIBias(i,:) = getConfInterval(biasDiffAll(:,i), pctCI);
end

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
    
    
% bestCIStim(i,:) = getConfInterval(reshape(stimDiff,numSearchIters,numStimStats), pctCI);
% bestCIHist(i,:) = getConfInterval(histDiffRand(:,popSize == i), pctCI);
% bestCIBias(i,:) = getConfInterval(histDiffRand(:,popSize == i), pctCI);
    
end

for i = 2:numSearchIters
    [stimPval(i), signifValsStim1(i)] = ranksum(stimDiffAll(:,i), stimDiffRand(:,popSize == i), 'alpha', .01);
    [histPval(i), signifValsHist1(i)] = ranksum(histDiffAll(:,i),histDiffRand(:,popSize == i), 'alpha', .01);
    [biasPval(i), signifValsBias1(i)] = ranksum(biasDiffAll(:,i), biasDiffRand(:,popSize == i), 'alpha', .01);
    [stimPval(i), signifValsStim2(i)] = ranksum(stimDiffAll(:,i), stimDiffRand(:,popSize == i), 'alpha', .05);
    [histPval(i), signifValsHist2(i)] = ranksum(histDiffAll(:,i),histDiffRand(:,popSize == i), 'alpha', .05);
    [biasPval(i), signifValsBias2(i)] = ranksum(biasDiffAll(:,i), biasDiffRand(:,popSize == i), 'alpha', .05);
end
    
signifValsStim1 = find(signifValsStim1==1);
signifValsHist1 = find(signifValsHist1==1);
signifValsBias1 = find(signifValsBias1==1);

signifValsStim2 = find(signifValsStim2==1);
signifValsHist2 = find(signifValsHist2==1);
signifValsBias2 = find(signifValsBias2==1);

%% plot best diversity values vs randomly sampled ones
divPan = panel();
divPan.pack('h', 3);

validInds = 2:numSearchIters;
for i = 1:3
    divPan(i).select();
    hold on;
    if i == 1;
        monteCarloMean = meanDivStim(validInds)';
        monteCarloCI = meanCIStim(validInds,:);
        bestMean = bestAvgDivStim((validInds));
        bestStdErr = stdErrDivStim(validInds);
        bestCI = bestCIStim(validInds,:);
        bestSignifVals1 = signifValsStim1;
        axisVals = [1 numSearchIters+1 -.2 1.6];
        plotSignifVal1 = 1.5;
        bestSignifVals2 = signifValsStim2;
        plotSignifVal2 = 1.3;
    elseif i == 2;
        monteCarloMean = meanDivHist(validInds)';
        monteCarloCI = meanCIHist(validInds,:);
        bestMean = bestAvgDivHist((validInds));
        bestStdErr = stdErrDivHist(validInds);
        bestCI = bestCIHist(validInds,:);
        bestSignifVals1 = signifValsHist1;
        axisVals = [1 numSearchIters+1 -5 50];
        plotSignifVal1 = 45;
        bestSignifVals2 = signifValsHist2;
        plotSignifVal2 = 40;
    elseif i == 3
        monteCarloMean = meanDivBias(validInds)';
        monteCarloCI = meanCIBias(validInds,:);
        bestMean = bestAvgDivBias((validInds));
        bestStdErr = stdErrDivBias(validInds);
        bestCI = bestCIBias(validInds,:);
        bestSignifVals1 = signifValsBias1;
        axisVals = [1 numSearchIters+1 -3 20];
        plotSignifVal1 = [];
        bestSignifVals2 = signifValsBias2;
        plotSignifVal2 = [];
    end
    monteCarloErrVals = [monteCarloMean - monteCarloCI(:,1), monteCarloCI(:,2) - monteCarloMean];
    bestErrVals = [bestMean' - bestCI(:,1), bestCI(:,2) - bestMean'];
    
%     boundedline(validInds, monteCarloMean, monteCarloErrVals, 'g', validInds, bestMean, bestErrVals, 'b', 'transparency', .2);
    
    boundedline(validInds, monteCarloMean, monteCarloErrVals, 'g', validInds, bestMean, bestStdErr, 'b', 'transparency', .2);
%     boundedline(validInds, bestMean, bestStdErr, 'b', 'transparency', .2);
    plot(2:numSearchIters, zeros(1, numSearchIters-1), 'r'); % represents homo diversity
    plot(bestSignifVals1, plotSignifVal1, 'b*');
    plot(bestSignifVals2, plotSignifVal2, 'c.');
    plot(validInds, bestMean, 'b.', validInds, monteCarloMean, 'g.', 2:numSearchIters, zeros(1, numSearchIters-1), 'r.')
    axis(axisVals);
%     ylabel('Population diversity');
    set(gca, 'XTick', [2 10])
end
divPan.ylabel('Population diversity');
divPan.xlabel('Population size');
divPan.fontsize = 6;

%% count how many unique neurons per best pop there are, and compare to chance
for i = 1:numStimStats
    numUniqueCells(i) = numel(unique(bestPop(:,i)));
end
numUniqueCells = numUniqueCells(includeStims);
numStimStats = length(includeStims);

histEdges = 0:10;
[uniCnts, uniCenters] = histc(numUniqueCells,histEdges);

% try to get error bars on expected number via random sampling
tempHist = zeros(10+1, 20000);
for i = 1:20000
    for j = 1:numStimStats
        tempBestPop = randsample(length(cellModels), 10, 1);
        randNumUnique(j) = numel(unique(tempBestPop));
    end
    tempHist(:,i) = histc(randNumUnique, histEdges);
end

expHist = median(tempHist');
pctCI = 50;

for i = 1:length(histEdges)
    expCI(:,i) = getConfInterval(tempHist(i,:), pctCI);
end

signifVals = histEdges(uniCnts < expCI(1,:) | uniCnts > expCI(2,:));
meanSignifVal = ranksum(numUniqueCells, randNumUnique); % p = 1.6*10-4
% plot(histEdges, expHist)

monteCarloErrVals = [expHist - expCI(1,:); expCI(2,:) - expHist];

numUniquePan = panel();
numUniquePan.select();
hold on;
bar(histEdges, uniCnts)
boundedline(histEdges', expHist', monteCarloErrVals', 'g', 'transparency', .2);
axis([0 11 0 7]);
xlabel('Number unique neurons')
ylabel('Population count');
numUniquePan.fontsize = 6;

%% put figures together
close all
clf
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0 0 18.3 8.0]);

greedySearchFig = panel();

greedySearchFig.pack('v', [66 -1]);
row1 = greedySearchFig(1);
row2 = greedySearchFig(2);

row1.pack('h', [33 -1]);
neuronDistnPan = row1(2);

row2.pack('h', [75 -1]);
divPan = row2(1);
numUniquePan = row2(2);
divPan
numUniquePan

%% plot an example figure showing error for avg het and hom populations compared to best pop
cd C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\greedySearch\finSims\randComparison
load 'dataProcessed';
cd C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\greedySearch\finSims\sim2\
load('greedySearch12', 'bestPopValAll');
temp = bestPopValAll;
load('greedySearch11', 'bestPopValAll');
bestPopValAll(2,:) = temp;

figure;
reconErrorPan = panel();
reconErrorPan.pack('h',2);

for i = 1:2
    reconErrorPan(i).select();
    hold on;
    boundedline(1:numSearchIters , avgErr(:,i), avgStdErr(:,i), 'g', 'transparency', .2);
    plot(1:numSearchIters, bestPopValAll(i,:), '-b');
end
reconErrorPan.ylabel('Decoding error (rmse)');
reconErrorPan.xlabel('Population size');
reconErrorPan.fontsize = 6;