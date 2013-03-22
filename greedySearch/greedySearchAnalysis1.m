%% analyze output of greedySearch

% figure out neuron addition path
currBestPop = [];
addedNeuron = [];
for i = 1:numSearchIters
    addedNeuron(i) = multisetdiff(bestPopListAll(i).dat, currBestPop);
    currBestPop = bestPopListAll(i).dat;
end
finBestPop = currBestPop;
for i = 1:numSearchIters
    numNeuronsInBest(i) = sum(addedNeuron(i)==finBestPop);
end
%% make plots showing addition of marginal neuron in glm parm space


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

data = kMat;
kMatMean = mean(kMat,2);
pcaEigs = 2;
[outputPCsK,pMatrixK, rankEigen] = princomp(data');
pctVarExpK = cumsum(rankEigen)./sum(rankEigen);

data = ihMat;
ihMatMean = mean(ihMat,2);
pcaEigs = 2;
[outputPCsH,pMatrixH, rankEigen] = princomp(data');
pctVarExpH = cumsum(rankEigen)./sum(rankEigen);

data = dcComp;
dcCompMean = mean(dcComp);
pcaEigs = 1;
[outputPCsB,pMatrixB, rankEigen] = princomp(data');
pctVarExpB = cumsum(rankEigen)./sum(rankEigen);

colorCode = colormap(jet(length(keepCells)));
pMatrix = pMatrixK;
outputPCs = outputPCsK;
pcaEigs = 2;
% figure; plot(outputPCs(1:50,1:pcaEigs));
figure; 
subplot(1,5,1:2);
hold on;
scatter(pMatrix(:,1), pMatrix(:,2), 12, colorCode, 'filled');
for i = 1:numSearchIters
    currNeuron = addedNeuron(i);
    text(pMatrix(currNeuron,1)+.05,pMatrix(currNeuron,2), num2str(i));
end

subplot(1,5,3:4);
hold on;
pMatrix = pMatrixH;
outputPCs = outputPCsH;
pcaEigs = 2;
scatter(pMatrix(:,1), pMatrix(:,2), 12, colorCode, 'filled');
for i = 1:numSearchIters
    currNeuron = addedNeuron(i);
    text(pMatrix(currNeuron,1)+2,pMatrix(currNeuron,2), num2str(i));
end

subplot(1,5,5);
hold on;
pMatrix = pMatrixB;
outputPCs = outputPCsB;
pcaEigs = 1;
scatter(ones(numModels,1), pMatrix(:,1), 12, colorCode, 'filled');
for i = 1:numSearchIters
    currNeuron = addedNeuron(i);
    text(ones(1,1)+.2,pMatrix(currNeuron,1), num2str(i));
end

% 
% figure; plot(tptsK, outputPCs(:,1), tptsK, outputPCs(:,2));

%% 

allFilts = [kMat' downsample(ihMat,10)' dcComp']';
allFilts = zscore(allFilts);

data = allFilts;
allFiltsMean = mean(allFilts,2);
pcaEigs = 2;
[outputPCsA,pMatrixA, rankEigen] = princomp(data');
pctVarExpA = cumsum(rankEigen)./sum(rankEigen);

pMatrix = pMatrixA;
outputPCs = outputPCsA;
colorCode = colormap(jet(length(keepCells)));
figure;
scatter(pMatrix(:,1), pMatrix(:,2), 12, colorCode, 'filled');
for i = 1:numSearchIters
    currNeuron = addedNeuron(i);
    text(pMatrix(currNeuron,1)+.5, pMatrix(currNeuron,2), num2str(i));
end
xlabel('PCA 1');
ylabel('PCA 2');

figure; hold on;
for i = 1:size(pMatrix,1)
    plot(pMatrix(i,1), pMatrix(i,2), '.', 'Color', colorCode(i,:,:));
end

for i = 1:numSearchIters
    currNeuron = addedNeuron(i);
    text(pMatrix(currNeuron,1)+.5, pMatrix(currNeuron,2), num2str(i));
end
xlabel('PCA 1');
ylabel('PCA 2');

% figure; plot(1:size(allFilts,1), outputPCs(:,1), 1:size(allFilts,1), outputPCs(:,2));
%% plot an example figure showing error for avg het and hom populations compared to best pop
cd 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\diffCorrs\diffDiversity\diffDiversity2\sim7'
load 'corrExptResults';
figure; hold on;
errorbar(popSizes(1:6), meanErrorsPerDivHet(1:6, 10), stdErrorsPerDivHet(1:6, 10), 'g');
errorbar(popSizes(1:6), meanErrorsPerDivHom(1:6, 10), stdErrorsPerDivHom(1:6, 10), 'r');
plot(1:numSearchIters, bestPopValAll, '.b')
ylabel('Recon error (rmse)');
xlabel('Population size');
axis([0 11 .4 .9])
set(gca, 'XTick', [1 10])