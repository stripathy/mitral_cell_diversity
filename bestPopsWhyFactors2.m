
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

% compute diversity among parameters
for i = 1:numPops
    [stimDiff(i) histDiff(i) biasDiff(i)] = glmDiff(cellModels(popList(:,i)));
end
diffMat = [stimDiff' histDiff' biasDiff'];


load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\drugAllFitsFinalggAll.mat';

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

allFilts = [kMat' downsample(ihMat,10)' dcComp']';
allFiltsMean = mean(allFilts,2);
allFiltsStd = std(allFilts')';

allFilts = zscore(allFilts);

data = allFilts;
% allFiltsMean = mean(allFilts,2);
pcaEigs = 2;
[outputPCsA,pMatrixA, rankEigen] = princomp(data');
pctVarExpA = cumsum(rankEigen)./sum(rankEigen);

pMatrix = pMatrixA;
outputPCs = outputPCsA;
colorCode = colormap(jet(length(keepCells)));
figure;
scatter(pMatrix(:,1), pMatrix(:,2), 12, colorCode, 'filled');
% scatter3(pMatrix(:,1), pMatrix(:,2), pMatrix(:,3), 12, colorCode, 'filled');

xlabel('PCA 1');
ylabel('PCA 2');

%% for each population, figure out it's contribution along each PCA


projFilters = [pMatrixA(:,1:3)];

projFilterMeans = mean(projFilters);
projFilterStds = std(projFilters);

projFilters = zscore(projFilters);

clear *FiltParms
for i = 1:numPops
    for j = 1:maxPopSize
        filtParmInds(j) = find(popList(j,i) == neuronsToSim);
    end
    meanProjFiltParms(i,:) = mean(projFilters(filtParmInds,:),1);
    stdProjFiltParms(i,:) = std(projFilters(filtParmInds,:),1);
    if numel(unique(filtParmInds)) > 1
        skewProjFiltParms(i,:) = skewness(projFilters(filtParmInds,:),[], 1);
    else
        skewProjFiltParms(i,:) = nan(1, size(projFilters,2));
    end
end

% Y = popAvgVal(:,1);
% X = [ones(numPops, 1) meanProjFiltParms];
% [b,bint,r,rint,stats] = regress(Y,X);
% rsqVals = stats(1)

% Y = zscore(mean(popAvgVal(:,1:10)')');

meanMeanParms = mean(meanProjFiltParms);
stdMeanParms = std(meanProjFiltParms);
meanStdParms = mean(stdProjFiltParms);
stdStdParms = std(stdProjFiltParms);

Y = (popAvgVal(:,1));
X = [ones(numPops, 1) (meanProjFiltParms) (stdProjFiltParms)];

Y = zscore(popAvgVal(:,2));
X = [ones(numPops, 1) zscore(meanProjFiltParms) zscore(stdProjFiltParms)];
[b,bint,r,rint,stats] = regress(Y,X);
rsqVals = stats(1)

bar(b);

testStims = [1 2 4 10];
for i = 1:length(testStims)
    subplot(2, 2, i);
    Y = zscore(popAvgVal(:,testStims(i)));
    X = [ones(numPops, 1) zscore(meanProjFiltParms) zscore(stdProjFiltParms)];
    [b,bint,r,rint,stats] = regress(Y,X);
    
    bar(b);
    betaCoeffs(i,:) = b;
    axis([0 size(X,2)+1 -.8 .5]);
end
% rsqVals = stats(1)

%% try converting standardized beta coeffs back into glm parms

i = 1;
meanGLMParms1 = ((meanMeanParms + betaCoeffs(i,2:4).*stdMeanParms).*projFilterStds)*outputPCsA(:,1:3)' + allFiltsMean';

i = 4;
meanGLMParms2 = ((meanMeanParms + betaCoeffs(i,2:4).*stdMeanParms).*projFilterStds)*outputPCsA(:,1:3)' + allFiltsMean';

for i = 1:length(testStims)
    meanBCoeffs = -betaCoeffs(i,2:4);
    a = stdMeanParms.*meanBCoeffs + meanMeanParms;
    b = a.*projFilterStds + projFilterMeans;
    c = b*outputPCsA(:,1:3)';
    d = c.*allFiltsStd' + allFiltsMean';
    bestFilts(i,:) = d;
end

d = c + allFiltsMean';

%% 
i = 3;
meanBCoeffs = -betaCoeffs(i,2:4);
a = stdMeanParms.*meanBCoeffs + meanMeanParms;

