%% look at data from the big pop sims, ask whether a certain population is
% the best

%% load data
cd 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\bigPopSims1\run_4\';
load 'dataPostAnalysis1';

%% find temp errors for each of the relevant pops

currPopSize = 5;
%find indexes of populations of currPopSize

cnt = 0;
for i = 1:length(pop)
    if length(pop(i).dat) == currPopSize
        relPopInds(cnt+1) = i;
        cnt = cnt + 1;
    end
end

popList = zeros(currPopSize, length(relPopInds));
for i = 1:length(relPopInds)
    popList(:,i) = pop(relPopInds(i)).dat;
end
%for each of the pops, compute their rmse's for each of the stims
stimEnt = .5*sum(myLog2(eig(stimCovMat))) + .5*slen*myLog2(2*pi*exp(1));
for i = 1:length(relPopInds)
%     hessDetSample = zeros(1, numSignals);
    for j = 1:numSignals
        meanReconAccSample(j) = sqrt(mean((outSignal(:,j) - reconStructOut(j+(relPopInds(i)-1)*numSignals).optStim).^2));
    end
    reconRMSEAll(:,i) = meanReconAccSample;
    reconRMSEMean(i) = mean(meanReconAccSample);
    
    respEntSample = zeros(1, numSignals);
    for j = 1:numSignals
%         hessMatTemp = reconStructOut(j+(i-1)*numSignals).hessian;
        respEntSample(j) = -.5*sum(myLog2(reconStructOut(j+(relPopInds(i)-1)*numSignals).hessianEigs));
    end

    respEnt(i) = mean(respEntSample) + .5*slen*myLog2(2*pi*exp(1));
    respEntError(i) = std(respEntSample)/sqrt(numSignals);
    totalEnt(i) = stimEnt - respEnt(i);
    totalEntError(i) = respEntError(i);
end

% which pop has the lowest mean error?
[sortedRMSEVals, bestPops] = sort(reconRMSEMean);

% do the same pops have consistently the lowest error?
[sortedRMSEAllVals, sortedBestCellIndices] = sort(reconRMSEAll, 2);

imagesc(sortedBestCellIndices);

popList(:,bestPops(1:10))

[sorted, sortedIndices] = sort(reconRMSEAll, 2);

[sortedEntVals, bestPopsByEnt] = sort(totalEnt, 'descend');
% reconRMSEOrder = 

%% find temp errors for each of the relevant pops when multiple stims are presented

currPopSize = 5;
%find indexes of populations of currPopSize

cnt = 0;
for i = 1:length(pop)
    if length(pop(i).dat) == currPopSize
        relPopInds(cnt+1) = i;
        cnt = cnt + 1;
    end
end

popList = zeros(currPopSize, length(relPopInds));
for i = 1:length(relPopInds)
    popList(:,i) = pop(relPopInds(i)).dat;
end
%for each of the pops, compute their rmse's for each of the stims
% stimEnt = .5*sum(myLog2(eig(stimCovMat))) + .5*slen*myLog2(2*pi*exp(1));
for i = 1:length(relPopInds)
%     hessDetSample = zeros(1, numSignals);
    for j = 1:numSignals
        meanReconAccSample(j) = sqrt(mean((outSignal(:,j) - reconStructOut(j+(relPopInds(i)-1)*numSignals).optStim).^2));
    end
    reconRMSEAll(:,i) = meanReconAccSample;
%     reconRMSEMean(i) = mean(meanReconAccSample);
%     
%     respEntSample = zeros(1, numSignals);
%     for j = 1:numSignals
% %         hessMatTemp = reconStructOut(j+(i-1)*numSignals).hessian;
%         respEntSample(j) = -.5*sum(myLog2(reconStructOut(j+(relPopInds(i)-1)*numSignals).hessianEigs));
%     end
% 
%     respEnt(i) = mean(respEntSample) + .5*slen*myLog2(2*pi*exp(1));
%     respEntError(i) = std(respEntSample)/sqrt(numSignals);
%     totalEnt(i) = stimEnt - respEnt(i);
%     totalEntError(i) = respEntError(i);
end

for i = 1:numStimStats
    reconRMSEMean(i,:) = mean(reconRMSEAll(1+(i-1)*numSignalsPerStat:i*numSignalsPerStat,:));
end
[sortedRMSEVals, bestPops] = sort(reconRMSEMean,2);


    figure;
    hist(nonzeros(popList(:,bestPops(1,1:10))), 11);
a = popList(:,bestPops(3,1:20));
for i = 1:3
    figure;
    hist(popList(:,sortedRMSEVals(1,1:10)));
end

% which pop has the lowest mean error?
[sortedRMSEVals, bestPops] = sort(reconRMSEMean);

% do the same pops have consistently the lowest error?
[sortedRMSEAllVals, sortedBestCellIndices] = sort(reconRMSEAll, 2);

imagesc(sortedBestCellIndices);

popList(:,bestPops(1:10))

[sorted, sortedIndices] = sort(reconRMSEAll, 2);

[sortedEntVals, bestPopsByEnt] = sort(totalEnt, 'descend');
% reconRMSEOrder = 
figure;
hist(nonzeros(bestPops(:,1:10)),numPops)
xlabel('population index'); ylabel('number occurences in top 10 pops');

% figure out how likely it is to get 4 pops which are best across 3 stims
popIndexThresh = 5; % the threshold to look for "good pops"
for i = 1:10000
    for j = 1:numStimStats
        tempPerm(j,:) = randperm(numPops);
    end
    [histCounts] = histc(nonzeros(tempPerm(:,1:popIndexThresh)), 0:1:numPops);
    numberTwos(i) = sum(histCounts==2);
    numberThrees(i) = sum(histCounts==3);
end

%% what is it about the best populations across the stims?

popIndexThresh = 10; % the threshold to look for "good pops"

%find the bestest pops
bestPopsAllStims = find(histc(nonzeros(bestPops(:,1:popIndexThresh)),0:1:numPops+1) >= 4)-1;

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
[count, bins] = histc(nonzeros(bestPopList(:,1:2)),0:1:length(neuronsToSim));
bar(bins, count);
figure; 
xlabel('neuron index'); ylabel('counts');

% plot the glm parms for a few pops
tempPopInd = 1;
% currCellInds = bestPopList(:,tempPopInd);
currCellInds = popList(:,624);
% currCellInds = [4 30];
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

%% plot bar graph showing which pops are best
threshPct = .05;
bestPopsThresh = floor(numPops*threshPct);

figure; 
for i = 1:numStimStats
    subplot(numStimStats+1,1,i);
    numHistBins = numPops;
    bar_handle = bar(bestPops(i,:), 1:numPops, 'BaseValue', numPops+1 , 'FaceColor', 'k');
    set(gca,'YDir','reverse');
    axis([0 numPops+1 0 numPops]);
    hold on;
    bar(bestPops(i,:), [1:bestPopsThresh ones(1, numPops-bestPopsThresh)*(numPops+2)], 'FaceColor', 'r');
    set(gca,'XTick', []);
end

%sort each bestest pop in order of lowest sum position across each stim
tempSum= zeros(numPops, 1);
for j = 1:numPops
    for i = 1:numStimStats
        tempSum(j) = tempSum(j) + find(bestPops(i,:) == j);
%         a(i,j) = find(bestPops(i,:) == j);
    end
end
meanPosition = tempSum/numStimStats;
subplot(numStimStats+1,1,numStimStats+1);
bar_handle = bar(1:numPops, meanPosition, 'BaseValue', max(meanPosition), 'FaceColor', 'k');
set(gca,'YDir','reverse');

%find the bestest pops
for i = 1:numStimStats
    tempVec = find(histc(nonzeros(bestPops(:,1:bestPopsThresh)),0:1:numPops+1) == i)-1;
    bestPopsAllStims(i, 1:length(tempVec)) = tempVec;
end


[numTimesInBest] = hist(nonzeros(bestPops(:,1:bestPopsThresh)),1:numPops);

colorOrder = ['r', 'b', 'g'];
hold on;

for i = 1:numStimStats
    tempVec = find(numTimesInBest==i);
    tempVecBar = ones(1,numPops)*(max(meanPosition)+2);
    tempVecBar(tempVec) = meanPosition(tempVec);
    bar(1:numPops, tempVecBar, 'FaceColor', colorOrder(i));
end
axis([0 numPops+1 0 max(meanPosition)]);
xlabel('Population index');
ylabel('Ranking');

%% 
sortedPopPositions = zeros(size(bestPops));
for i = 1:numStimStats
    for j = 1:numPops
        sortedPopPositions(i,bestPops(i,j)) = j;
    end
end

[sortA sortedAInds] = sort(sortedPopPositions(1,:));

figure;
plot(sortA, sortedPopPositions(3,sortedAInds),'.')
set(gca,'YDir','reverse','XDir','reverse');
axis ([-10 904 -10 904]);




% diff between freq selectivity of neurons in pop?


    
    




