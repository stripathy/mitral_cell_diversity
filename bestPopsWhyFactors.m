% try to explain why best pops are the best in terms of statistics
% computed on pops

%% first compute or use already computed population statistics
meanFRBoth = mean(meanFRPerStim,2);
corrMatBoth = mean(corrMat, 3);
reconRMSEMean(3,:) = mean(reconRMSEMean);

meanFRPerStim(:,3) = meanFRBoth;
corrMat(:,:,3) = corrMatBoth;

totalEntPerSec = totalEntPerSec';
totalEntPerSec(:,3) = mean(totalEntPerSec,2);

% compute mean FR, corr, reliab
for k = 1:3
    stimInd = k;
    for i = 1:numPops
        popMeanFR(i,k) = mean(meanFRPerStim(popList(:,i),stimInd));
        popMeanInfo(i,k) = mean(totalEntPerSec(popList(:,i),stimInd));
        tempPerms = combnk(popList(:,i),2);
        for j = 1:length(tempPerms)
            popMeanCorrAll(j,i) = corrMat(tempPerms(j,1), tempPerms(j,2), stimInd);
        end
    end
    popMeanCorr(:,k) = mean(popMeanCorrAll);
    
    meanReliab = diag(corrMat(:,:, stimInd));
    for i = 1:numPops
        popMeanReliab(i,k) = mean(meanReliab(popList(:,i)));
    end
end

% compute diversity among parameters
for i = 1:numPops
    [stimDiff(i) histDiff(i) biasDiff(i)] = glmDiff(cellModels(popList(:,i)));
end
diffMat = [stimDiff' histDiff' biasDiff'];


% numbers of neurons in each population
popMakeup = histc(popList,1:44);

dataMat = [reconRMSEMean(1,:)' reconRMSEMean(2,:)' popMeanFR' popMeanInfo' popMeanCorr' popMeanReliab' stimDiff' histDiff' biasDiff'];

dataLabels = {'errors1', 'errors2', 'FR', 'corr', 'reliab', 'stimDiff', 'histDiff', 'biasDiff'};

a = dataLabels;

X = [ones(numPops, 1) popMeanFR(:,3)];
Y = [reconRMSEMean(3,:)' reconRMSEMean(2,:)'];

X = [ones(numPops, 1) popMeanFR(:,3) stimDiff' histDiff' biasDiff'];

X11 = [ones(numPops, 1) popMeanFR(:,1)];
X12 = [ones(numPops, 1) popMeanFR(:,1) stimDiff' histDiff' biasDiff'];

frMat = ones(numPops, 1) popMeanFR(:,3)



for i = 1:3
    Y = reconRMSEMean(i,:)';
    for j = 1:2
        if j ==1
            X = [ones(numPops, 1) popMeanFR(:,i)];
            [b,bint,r,rint,stats] = regress(Y,X);
            rsqVals(i,j) = stats(1);
        else
            X = [ones(numPops, 1) popMeanFR(:,i) diffMat];
            [b,bint,r,rint,stats] = regress(Y,X);
            rsqVals(i,j) = stats(1);
        end
    end
end

%% compute pca projections for populations

histInds = 101:size(cellModels(1).iht,1);
for j = 1:length(neuronsToSim)
    currCell = neuronsToSim(j);
    kMat(:,j) = cellModels(currCell).k;
    histTemp = exp(cellModels(currCell).ihbas*cellModels(currCell).ih);
    ihMat(:,j) = log10(histTemp(histInds));
    dcComp(j) = exp(cellModels(currCell).dc);
end


allFilts = [kMat' ihMat' dcComp']';

tptsK = -50+1:1:0;

figure;
for j = 1:length(keepCells)
%             subplot(1,4,1); hold on;
%             plot(tptsK, staMat(:,j), 'Color', colorCode(j,:), 'LineWidth' , 2); axis([-50 0 -.5 1.5]); 
        subplot(1,3,1); hold on; plot(tptsK, kMat(:,j), 'Color',colorCode(j,:), 'LineWidth' , 1.5); axis tight; %axis([-50 0 -.5 1.2]); 
        ylabel('log(cond. intensity)'); xlabel('time (ms)');
        subplot(1,3,2); hold on; plot(ggAll(1).iht(histInds),ihMat(:,j), 'Color',colorCode(j,:), 'LineWidth' , 2); ylabel('Gain'); axis tight; 
        subplot(1,3,3); hold on; plot(1, dcComp(j), '.', 'Color',colorCode(j,:), 'LineWidth' , 2); set(gca, 'XTick', []); axis tight; 
end
legend(legendStr)
colorCode = colormap(jet(length(keepCells)));

numEigsK = 3;
numEigsH = 2;
numEigsB = 1;

data = kMat;
pcaEigs = 3;
[outputPCsK,pMatrixK, rankEigen] = princomp(data');
pctVarExpK = cumsum(rankEigen)./sum(rankEigen);

data = ihMat;
pcaEigs = 2;
[outputPCsH,pMatrixH, rankEigen] = princomp(data');
pctVarExpH = cumsum(rankEigen)./sum(rankEigen);

data = dcComp;
pcaEigs = 1;
[outputPCsB,pMatrixB, rankEigen] = princomp(data');
pctVarExpB = cumsum(rankEigen)./sum(rankEigen);

projFilters = [pMatrixK(:,1:numEigsK) pMatrixH(:,1:numEigsH) pMatrixB(:,1:numEigsB)];

projFilterMeans = mean(projFilters);
projFilterStds = std(projFilters);

projFilters = zscore(projFilters);


pMatrix = pMatrixH;
outputPCs = outputPCsH;
pcaEigs = numEigsH;
figure; plot(outputPCs(1:length(outputPCs),1:pcaEigs));
% [outputPCs, rankEigen, pMatrix] = kpPCA(kMat(:,keepCells), pcaEigs);
figure;  hold on;
for i = 1:size(pMatrix,1)
    plot(pMatrix(i,1), pMatrix(i,2),'*', 'Color',colorCode(i,:), 'LineWidth' , 2);
    text(pMatrix(i,1)+.05,pMatrix(i,2), num2str(keepCells(i)));
end

figure; plot(tptsK, outputPCs(:,1), tptsK, outputPCs(:,2));

% compute averages across projFilters for each pop
for i = 1:numPops
    for j = 1:maxPopSize
        filtParmInds(j) = find(popList(j,i) == neuronsToSim);
    end
    meanProjFiltParms(i,:) = mean(projFilters(filtParmInds,:),1);
    stdProjFiltParms(i,:) = std(projFilters(filtParmInds,:),1);
    skewProjFiltParms(i,:) = skewness(projFilters(filtParmInds,:),1);
end

Y = [zscore(reconRMSEMean(3,:)')];

X = [ones(numPops, 1) popMeanFR(:,1)];

X = [ones(numPops, 1) zscore(popMeanFR(:,3)) zscore(meanProjFiltParms) zscore(stdProjFiltParms)];
X = [ones(numPops, 1) popMeanFR(:,2) meanProjFiltParms stdProjFiltParms meanProjFiltParms.*meanProjFiltParms];

for i = 1:3
    Y = reconRMSEMean(i,:)';
    for j = 1:3
        if j ==1
            X = [ones(numPops, 1) popMeanInfo(:,i)];
            [b,bint,r,rint,stats] = regress(Y,X);
            rsqVals(i,j) = stats(1);
        elseif j == 2
            X = [ones(numPops, 1)  popMeanInfo(:,i) meanProjFiltParms];
            [b,bint,r,rint,stats] = regress(Y,X);
            rsqVals(i,j) = stats(1);
        elseif j == 3
            X = [ones(numPops, 1)  zscore([popMeanInfo(:,i) meanProjFiltParms stdProjFiltParms])];
            [b,bint,r,rint,stats] = regress(Y,X);
            rsqVals(i,j) = stats(1);
        elseif j == 4
            X = [ones(numPops, 1)  meanProjFiltParms stdProjFiltParms skewProjFiltParms];
            [b,bint,r,rint,stats] = regress(Y,X);
            rsqVals(i,j) = stats(1);
        elseif j == 5
            X = [ones(numPops, 1)  popMeanInfo(:,i) diffMat];
            [b,bint,r,rint,stats] = regress(Y,X);
            rsqVals(i,j) = stats(1);
        end
    end
end

dataMat = zscore([reconRMSEMean(1:3,:)' popMeanFR popMeanInfo meanProjFiltParms stdProjFiltParms]);

labels = y1 y2 y3 fr1 fr2 fr3 in1 in2 in3 spm1 spm2 spm3 hpm1 hpm2 bpm1 sps1 sps2 sps3 hps1 hps2 bps1

dataMat = zscore([reconRMSEMean(3,:)' popMeanFR(:,3) popMeanInfo(:,3) meanProjFiltParms stdProjFiltParms diffMat]);

labels = y3 fr3 in3 spm1 spm2 spm3 hpm1 hpm2 bpm1 sps1 sps2 sps3 hps1 hps2 bps1 sd hd bd

%% plot some stuff

% first plot a bar graph showing pct var explained for different models
% rsqValsEst = [77.6 81.5 89.9];
rsqValsEst = [63.1 67.6 89.8];

bar(rsqValsEst)
ylabel('Pct variance explained (R^2)');
axis([0 length(rsqValsEst)+1 50 100])

tptsK = -50+1:1:0;
pMatrix = pMatrixK;
outputPCs = outputPCsK;
pcaEigs = numEigsK;
figure; plot(tptsK, outputPCs(1:length(outputPCs),1:pcaEigs));
legend('stim pca1', 'stim pca2', 'stim pca3');
xlabel('Time (ms)');
ylabel('Stim value (a.u.)');
axis tight

pMatrix = pMatrixH;
outputPCs = outputPCsH;
pcaEigs = numEigsH;
figure; plot(ggAll(1).iht(histInds), outputPCs(1:length(outputPCs),1:pcaEigs));
xlabel('Time (ms)');
ylabel('Hist value (log gain)');
legend('hist pca1', 'hist pca2');
axis tight

%% find regression coeffs and conf ints

Yreal = Y;
Xreal = X;

numResamps = 1000;

for i = 1:numResamps
    sampInds = randi(numPops,numPops,1);
Xresamp = [ones(numPops, 1)  zscore([popMeanInfo(sampInds,3) meanProjFiltParms(sampInds,:) stdProjFiltParms(sampInds,:)])];

Yresamp = Yreal(sampInds,:);
b(i,:) = regress(Yresamp,Xresamp);
end

for i = 1:size(b,2)
    bConfInt(i,:) = getConfInterval(b(:,i));
    bMeans(i) = mean(b(:,i));
end

figure;
bar(bMeans)
hold on;
errorbar(1:length(bMeans), bMeans, bMeans' - bConfInt(:,1), bMeans' - bConfInt(:,2), '.k')

axis([1.5 14.5 -.8 .6])
xlabel(''); ylabel('Standardized regression coeff');
set(gca, 'XTick', []);

%% modify old data structures to only have high and low freq stim

totPops = 894;
popInds = 45:totPops;
popList = popList(:,popInds);
reconRMSEMean = reconRMSEMean([1 3],popInds);
numPops = length(popInds);
