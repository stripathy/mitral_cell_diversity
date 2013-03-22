
% first compute pca space of neurons
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\drugAllFitsFinalggAll.mat';
clear new* output* pM* pct* *Mat dc* pca* rank* *Mean

cellModels = ggAll;
keepCells = 1:44;
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

allFilts = [kMat' downsample(ihMat,10)' dcComp'];
allFiltsMean = mean(allFilts);
allFiltsStd = std(allFilts);

allFilts = zscore(allFilts);
% 
% allFilts = [kMat' downsample(ihMat,10)' dcComp']';
% allFiltsMean = mean(allFilts,2);
% allFiltsStd = std(allFilts,[],2);
% allFilts = zscore(allFilts);

data = allFilts;
% allFiltsMean = mean(allFilts,2);
pcaEigs = 2;
[outputPCsA,pMatrixA, rankEigen] = princomp(data);
pctVarExpA = cumsum(rankEigen)./sum(rankEigen);


pMatrix = pMatrixA;
outputPCs = outputPCsA;
colorCode = colormap(jet(length(keepCells)));

figure;
hold on;
for i = 1:size(pMatrix,1)
    plot(pMatrix(i,1), pMatrix(i,2), '*', 'Color', colorCode(i,:,:), 'MarkerSize', 10);
end
% scatter(pMatrix(:,1), pMatrix(:,2), sizeVec, colorCode, 'filled');

xlabel('PCA 1');
ylabel('PCA 2');

%%
clear new*  pM* pct* *Mat dc* pca* rank*
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\drugAllFitsFinalggAll.mat';

cellModels = ggAll;
keepCells = 1:44;
numModels = length(ggAll);

clear kMat ihMat
histInds = 1:length(ggAll(1).ihbas);
for j = 1:length(keepCells)
    currCellIds(j) = keepCells(j);%find(allGroupInds(i,j)==cellId);
%     staMat(:,j) = ggAll(currCellIds(j)).sta;
    kMat(:,j) = ggAll(currCellIds(j)).k;
    temp = log(exp(ggAll(currCellIds(j)).ihbas*ggAll(currCellIds(j)).ih));
    ihMat(:,j) = temp(histInds);
%     ihMat(:,j) = exp(ggAll(currCellIds(j)).ihbas*ggAll(currCellIds(j)).ih);
    dcComp(j) = log(exp(ggAll(currCellIds(j)).dc));
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

% divFact = [0 .5 .8 .9 .95 1 1.05 1.1 1.15 1.2 1.25 1.3];
% 
% divFact = [1];
% divFact = [1];
useComps = 20;
clear modCellModels
modelInd = 1;
divFact =.25;
numReps = 50;

cnt = 1;
for i = 1:numModels
    for j = 1:numReps
        stimDiv = (1 + randn*divFact);
        histDiv = (1 + randn*divFact);
        biasDiv = (1 + randn*divFact);
        if j == 1
            stimDiv = 1;
            histDiv = 1;
            biasDiv = 1;
        end
        
%         newKMat = repmat(kMatMean', numModels, 1) + pMatrixK(:,1:useComps)*outputPCsK(:,1:useComps)'*stimDiv;
%         newihMat = repmat(ihMatMean', numModels, 1) + pMatrixH(:,1:useComps)*outputPCsH(:,1:useComps)'*histDiv;
        newKMat = repmat(kMatMean', numModels, 1) + pMatrixK*outputPCsK'*stimDiv;
        newihMat = repmat(ihMatMean', numModels, 1) + pMatrixH*outputPCsH'*histDiv;
        newbMat = repmat(dcCompMean', numModels, 1) + pMatrixB*outputPCsB'*biasDiv;
        
        modCellModels(cnt) = ggAll(i);
        modCellModels(cnt).k = newKMat(i,:)';
        modCellModels(cnt).kt = ggAll(i).ktbas\newKMat(i,:)';
        
        modCellModels(cnt).ih = ggAll(i).ihbas\newihMat(i,:)';
        modCellModels(cnt).dc = newbMat(i);
        %                 modCellModels(modelInd).ih = ggAll(25).ih;
        %         modCellModels(modelInd).dc = ggAll(25).dc;
        modelInd = modelInd + 1;
        cnt = cnt + 1;
    end
end

clear new* pM* pct* *Mat dc* pca* rank*

%% 
% 
cellModels = ggAll;
keepCells = 1:44;
numModels = length(ggAll);

clear kMat ihMat
histInds = 101:length(ggAll(1).ihbas);
for j = 1:length(modCellModels)
    currCellIds(j) = j;%find(allGroupInds(i,j)==cellId);
%     staMat(:,j) = ggAll(currCellIds(j)).sta;
    kMat(:,j) = modCellModels(currCellIds(j)).k;
    temp = log10(exp(modCellModels(currCellIds(j)).ihbas*modCellModels(currCellIds(j)).ih));
    ihMat(:,j) = temp(histInds);
%     ihMat(:,j) = exp(ggAll(currCellIds(j)).ihbas*ggAll(currCellIds(j)).ih);
    dcComp(j) = log10(exp(modCellModels(currCellIds(j)).dc));
end
allFiltsMod = [kMat' downsample(ihMat,10)' dcComp'];
allFiltsMod = (allFiltsMod - repmat(allFiltsMean, length(modCellModels), 1))./repmat(allFiltsStd, length(modCellModels),1);

projData = allFiltsMod*(outputPCsA);
% 
% allFiltsMean = mean(allFilts,2);
% pcaEigs = 2;
% [outputPCsA,pMatrixA, rankEigen] = princomp(data');
% pctVarExpA = cumsum(rankEigen)./sum(rankEigen);


pMatrix = projData;
% outputPCs = outputPCsA;
colorCode = colormap(jet(length(keepCells)));
cnt = 1;
for j = 1:numModels
    for i = 1:numReps
        colorCodeNew(cnt, :) = colorCode(j,:);
        cnt = cnt + 1;
    end
end

% figure;
hold on;
for i = 1:size(pMatrix,1)
    plot(pMatrix(i,1), pMatrix(i,2), '.', 'Color', colorCodeNew(i,:,:), 'MarkerSize', 10);
end
% scatter(pMatrix(:,1), pMatrix(:,2), sizeVec, colorCode, 'filled');

xlabel('PCA 1');
ylabel('PCA 2');
%%
clear new* pM* pct* *Mat dc* pca* rank* output* 

%%
cnt = 1;
for j = 1:44
    for i = 1:numReps
        modelNeuronInd(cnt) = j;
        cnt = cnt + 1;
    end
end

keepCellCnt = 5;
keepModelInds = [];
for i = 1:44
    inds = find(modelNeuronInd == i);
    validFRs = meanFR(inds(1)).*[1-divFact 1+divFact];
    validModels = find(meanFR(inds(2:end)) > validFRs(1) &  meanFR(inds(2:end)) < validFRs(2));
    a(i) = length(validModels);
    keepInds = randsample(inds(validModels), keepCellCnt);
    keepInds(1) = inds(1);
    keepModelInds = [keepModelInds keepInds];
end

newCellModels = ggAll(keepModelInds);
oldCellModels = ggAll;