% best pops fig generation

% plot frequencies of each stimulus distn
tapers=3; avg=1; fs = 1000;
pars = struct ('Fs', fs, ...
    'tapers', [tapers 2*tapers-1], ...
    'pad', 0, ...
    'trialave', 1);
for k = 1:numStimStats
    [S(:,k), f] = mtspectrumc(outSignal(:,1+(k-1)*numSignalsPerStat:k*numSignalsPerStat), pars);
end
figure; plot(f,S); ylabel('Power spectra'); xlabel('Frequency (Hz)');
axis([0 100 0 .045]);
%% find temp errors for each of the relevant pops when multiple stims are presented

load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\diffStimStats\sim1\reconsMerged.mat'; 

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
[blah, meanFRSortedPops] = sort(sum(popList,1));
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
    reconRMSEMean(i,:) = mean(reconRMSEAll(1+(i-1)*numSignalsPerStat:i*numSignalsPerStat,1:numPops));
end
[sortedRMSEVals, bestPops] = sort(reconRMSEMean,2);

figure; 
bestPopsThresh = floor(numPops*.1);
for i = 1:numStimStats
    subplot(4,1,i);
    numHistBins = numPops;
    bar_handle = bar(bestPops(i,:), 1:numPops, 'BaseValue', numPops+1 , 'FaceColor', 'k');
    set(gca,'YDir','reverse');
    axis([0 numPops+1 0 numPops]);
    hold on;
    bar(bestPops(i,:), [1:bestPopsThresh ones(1, numPops-bestPopsThresh)*(numPops+2)], 'FaceColor', 'r');
    set(gca,'XTick', []); box off;
    ylabel('Ranking');
end

%sort each bestest pop in order of lowest sum position across each stim
tempSum= zeros(numPops, 1);
for j = 1:numPops
    for i = 1:numStimStats
        tempSum(j) = tempSum(j) + find(bestPops(i,:) == j);
%         a(i,j) = find(bestPops(i,:) == j);
    end
end
meanPosition = tempSum/3;
subplot(4,1,4);
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
ylabel('Avg. Ranking');


% plot relative positions of stim 1 vs stim 3
sortedPopPositions = zeros(size(bestPops));
for i = 1:numStimStats
    for j = 1:numPops
        sortedPopPositions(i,bestPops(i,j)) = j;
    end
end

stimInd1 = 1;
stimInd2 = 3;
[sortA sortedAInds] = sort(sortedPopPositions(stimInd1,:));

figure;
hold on;
plot([1 numPops], [1 numPops], 'k');
plot(sortA(heteroStimAInds), heteroStimBInds,'k.')

plot(sortA(homoStimAInds), homoStimBInds, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', [.5 .5 .5], 'LineWidth', 1)

plot(sortA(homoPopAInd), homoPopBInd,'ro', 'MarkerSize', 4, 'LineWidth', 1)
plot(sortA(heteroPopAInd), heteroPopBInd,'g.', 'MarkerSize', 10)

plot(sortA(bestPopHighAInd), bestPopHighBInd,'rh', 'MarkerSize', 10)
plot(sortA(bestPopLowAInd), bestPopLowBInd,'bh', 'MarkerSize', 10)
plot(sortA(bestPopBothAInd), bestPopBothBInd,'gh', 'MarkerSize', 10)
% plot(sortA, sortedPopPositions(stimInd2,sortedAInds),'k.')

set(gca,'YDir','reverse','XDir','reverse');
axis ([-10 904 -10 904]);
xlabel('Population rank: stim 1');
ylabel('Population rank: stim 2');

% plot(1:44, homoInds, 'ro')

homoStimAInds = sortedPopPositions(stimInd1,1:43);
homoStimBInds = sortedPopPositions(stimInd2,1:43);

heteroStimAInds = sortedPopPositions(stimInd1,45:end);
heteroStimBInds = sortedPopPositions(stimInd2,45:end);

homoPopAInd = sortedPopPositions(stimInd1, 44);
homoPopBInd = sortedPopPositions(stimInd2, 44);
heteroPopAInd = sortedPopPositions(stimInd1, 624);
heteroPopBInd = sortedPopPositions(stimInd2, 624);

bestPopHighAInd = sortedPopPositions(stimInd1, 895);
bestPopHighBInd = sortedPopPositions(stimInd2, 895);
bestPopLowAInd = sortedPopPositions(stimInd1, 896);
bestPopLowBInd = sortedPopPositions(stimInd2, 896);
bestPopBothAInd = sortedPopPositions(stimInd1, 898);
bestPopBothBInd = sortedPopPositions(stimInd2, 898);

for i = 1:44
%     homoInds(i) = find(sortedAInds==i);
    homoStimAInds(i) = find(sortedPopPositions(stimInd1,:) == i);
    homoStimBInds(i) = find(sortedPopPositions(stimInd2,:) == i);
end

% how many specialists and generalists would you expect by chance?

% first count how many pops are in top X % empirically
bestPopsThresh = floor(.1*numPops);
spec1 = numel(setdiff(bestPops(stimInd1,1:bestPopsThresh),  bestPops(stimInd2,1:bestPopsThresh)));
spec2 = numel(setdiff(bestPops(stimInd2,1:bestPopsThresh),  bestPops(stimInd1,1:bestPopsThresh)));
general = numel(intersect(bestPops(stimInd2,1:bestPopsThresh),  bestPops(stimInd1,1:bestPopsThresh)));
neither = numPops - (spec1 + spec2 + general);

numResamples = 1000;
for i = 1:numResamples
    randInts1 = randperm(numPops);
    randInts2 = randperm(numPops);
    resampSpec1(i) = numel(setdiff(randInts1(1:bestPopsThresh), randInts2(1:bestPopsThresh)));
    resampSpec2(i) = numel(setdiff(randInts2(1:bestPopsThresh), randInts1(1:bestPopsThresh)));
    resampGen(i) = numel(intersect(randInts1(1:bestPopsThresh), randInts2(1:bestPopsThresh)));
    resampNeither(i) = numPops - resampSpec1(i) - resampSpec2(i) - resampGen(i);
end
pctCI = 5;
confIntResampSpec1 = getConfInterval(resampSpec1, pctCI);
confIntResampSpec2 = getConfInterval(resampSpec2, pctCI);
confIntResampGen = getConfInterval(resampGen, pctCI);
confIntResampNeither = getConfInterval(resampNeither, pctCI);
    
barMat = [neither spec1 spec2 general; ...
    mean(confIntResampNeither), mean(confIntResampSpec1), mean(confIntResampSpec2), mean(confIntResampGen)];
barMat = reshape(barMat, 8,1);

confIntMeans = [mean(confIntResampNeither), mean(confIntResampSpec1), mean(confIntResampSpec2), mean(confIntResampGen)];
confInts = [confIntResampNeither' confIntResampSpec1' confIntResampSpec2' confIntResampGen'];
% plot the bar graph
figure;
bar(barMat')
hold on;
errorbar([2:2:8], confIntMeans, confIntMeans - confInts(1,:), '.k')
%% plot some example recons between 
currStimInd = 156;%20;
popInds = [36 196];
% popInds = [44 624];
for i = 1:size(popInds,2)
    plotInds(i) = currStimInd + numSignals*(popInds(i)-1);
end

colorCodePlot(1,:) = [1 0 0];
colorCodePlot(2,:) = [0 1 0];
colorCode = colormap(jet(length(keepCells)));

figure; subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2) % axis([1 slen -3 3]);
testSpikesPlot = [];
hold on;
rasterColors = [];
for i = 1:2%length(currCellInds)
    tempCurrCellInds = reconStructIn(plotInds(i)).popMakeup;
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim, 'Color', colorCodePlot(i,:), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    
%     rasterColors(end+1:end+nnz(tempCurrCellInds),:) = repmat(colorCodePlot(i,:), nnz(tempCurrCellInds),1);
    rasterColors(end+1:end+nnz(tempCurrCellInds),:) = colorCode(tempCurrCellInds,:);
    if i ~= 2
        rasterColors(end+1,:) = [1,1,1]; %white
        testSpikesPlot = vectCat(testSpikesPlot, 0);
    end
end
axis([0 slen -3 3]);
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('Stimulus (a.u.)'); set(gca,'XTickLabel', []); axis tight; box off;
subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('Time (ms)'); set(gca,'YTick', []);

%% plot some example recons and rasters for some of the pops
currStimInd = 135;

popInds = [38 20];
for i = 1:size(popInds,2)
    plotInds(i) = currStimInd + numSignals*(popInds(i)-1);
end

colorCodePlot(1,:) = [0 1 0];
colorCodePlot(2,:) = [1 0 0];

figure; subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2) % axis([1 slen -3 3]);
testSpikesPlot = [];
hold on;
rasterColors = [];
for i = 1:2%length(currCellInds)
    tempCurrCellInds = reconStructIn(plotInds(i)).popMakeup;
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim, 'Color', colorCodePlot(i,:), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    
    rasterColors(end+1:end+nnz(tempCurrCellInds),:) = repmat(colorCodePlot(i,:), nnz(tempCurrCellInds),1);
%     rasterColors(end+1:end+nnz(tempCurrCellInds),:) = colorCode(tempCurrCellInds,:);
    if i ~= 2
        rasterColors(end+1,:) = [1,1,1]; %white
        testSpikesPlot = vectCat(testSpikesPlot, 0);
    end
end
axis([0 slen -3 3]);
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('Stimulus (a.u.)'); set(gca,'XTickLabel', []); axis tight; box off;
subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('Time (ms)'); set(gca,'YTick', []);

%% remake an example of the 
[includePops] = [44 624];
otherPopsHomo = randsample(1:43, 9);
otherPopsHetero = randsample(45:numPops, 39);

includePops = sort([includePops otherPopsHomo otherPopsHetero]);
for i = 1:length(includePops)
    for j = 1:numStimStats
        popPosSmall(j,i) = find(bestPops(j,:) == includePops(i));
    end
end

[vals bestPopsSmall] = sort(popPosSmall, 2);
% bestPopsSmall = bestPops(:,includePops);

firstNeuronLoc = find(includePops == 44);
secondNeuronLoc = find(includePops == 624);

for i = 1:numStimStats
    firstNeuronRank(i) = find(bestPopsSmall(i,:) == firstNeuronLoc);
    secondNeuronRank(i) = find(bestPopsSmall(i,:) == secondNeuronLoc);
end

numPopsSmall = length(includePops);

plotInds = [1 3];
figure; 
bestPopsThresh = floor(numPops*.1);
for i = [1 2]%1:numStimStats
    subplot(2,1,i);
    numHistBins = numPopsSmall;
    bar_handle = bar(bestPopsSmall(plotInds(i),:), 1:numPopsSmall, 'BaseValue', numPopsSmall+1 , 'FaceColor', 'k');
    set(gca,'YDir','reverse');
    axis([0 numPopsSmall 0 numPopsSmall]);
    hold on;
    bar(firstNeuronLoc, firstNeuronRank(plotInds(i)), 'FaceColor', 'r');
    bar(secondNeuronLoc, secondNeuronRank(plotInds(i)), 'FaceColor', 'g');
%     bar(bestPopsSmall(plotInds(i),bestPopsSmall(plotInds(i),:)==find(includePops==624)), find(bestPopsSmall(plotInds(i),:)==bestPopsSmall(plotInds(i),bestPopsSmall(plotInds(i),:)==find(includePops==624))), 'FaceColor', 'g') 
%     bestPopsSmall(plotInds(i),find(includePops==44)), 'FaceColor', 'r');
%     bar(bestPops(i,:), [1:bestPopsThresh ones(1, numPops-bestPopsThresh)*(numPops+2)], 'FaceColor', 'r');
    set(gca,'XTick', []); box off;
    ylabel(['Ranking: stim ' num2str(i)]);
end

%sort each bestest pop in order of lowest sum position across each stim
tempSum= zeros(numPops, 1);
for j = 1:numPops
    for i = 1:numStimStats
        tempSum(j) = tempSum(j) + find(bestPops(i,:) == j);
%         a(i,j) = find(bestPops(i,:) == j);
    end
end
meanPosition = tempSum/3;
subplot(4,1,4);
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
ylabel('Avg. Ranking');

%% plot the glm parms for a few pops
tempPopInd = 1;
% currCellInds = bestPopList(:,tempPopInd);
currCellInds = popList(:,44);
% currCellInds = [4 30];
colorCode = colormap(jet(length(keepCells)));

figure; 
box off;
allCellInds = unique(nonzeros(currCellInds));
for i = 1:length(allCellInds)
    
    subplot(1,5,1:2);plot(-49:0, cellModels(allCellInds(i)).k, 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 2); axis([-50 0 -.2 1]); hold on;
    subplot(1,5,3:4);plot(cellModels(1).iht,log10(exp(cellModels(allCellInds(i)).ihbas*cellModels(allCellInds(i)).ih)), 'Color', colorCode(allCellInds(i),:), 'LineWidth' , 2); hold on;
    axis ([0 60 -1.6 2.2]);
    yTickPos = [-1:2]; yTickLabels = [.1 1 10 100];
    set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabels);
    subplot(1,5,5);plot(1, log10(exp(cellModels(allCellInds(i)).dc)),  '*', 'Color', colorCode(allCellInds(i),:)); hold on;
    axis ([0 2 -1.6 2.2]);
    set(gca, 'YTick', yTickPos, 'YTickLabel', [], 'XTick', []);
end
subplot(1,5,1:2);box off;
subplot(1,5,3:4); box off;
subplot(1,5,5);box off;