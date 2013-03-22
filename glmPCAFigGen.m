%% for an example stim stat, plot the iterative ordering of the neurons

% first compute pca space of neurons
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
    temp = log10(exp(ggAll(currCellIds(j)).ihbas*ggAll(currCellIds(j)).ih));
    ihMat(:,j) = temp(histInds);
%     ihMat(:,j) = exp(ggAll(currCellIds(j)).ihbas*ggAll(currCellIds(j)).ih);
    dcComp(j) = log10(exp(ggAll(currCellIds(j)).dc));
end
%% plot the glm parms for a few pops
popInds = 155;
figure;
glmParmPan = panel();
glmParmPan.pack('v', 1);

for j = 1:length(popInds)
    currGlmPan = glmParmPan(j);
    currGlmPan.pack('h', [45 45 -1]);
    kFilt = currGlmPan(1);
    hFilt = currGlmPan(2);
    bFilt = currGlmPan(3);
    
currCellInds = pop(popInds(j)).dat;


% figure; 
% box off;
colorCode = colormap(jet(length(keepCells)));
allCellInds = unique(nonzeros(currCellInds));
for i = 1:length(allCellInds)
    
    kFilt.select(); hold on; plot(-49:0, cellModels(allCellInds(i)).k, 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 1); 
    axis([-50 0 -.2 1.1]); xTickPos = [-50 -25 0]; hold on;
    set(gca, 'XTick', xTickPos);
    hFilt.select(); plot(cellModels(1).iht,log10(exp(cellModels(allCellInds(i)).ihbas*cellModels(allCellInds(i)).ih)), 'Color', colorCode(allCellInds(i),:), 'LineWidth' , 1); hold on;
    axis ([0 60 -1.6 2.2]);
    yTickPos = [-1:2]; yTickLabels = [.1 1 10 100];
    set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabels);
    bFilt.select(); plot(1, log10(exp(cellModels(allCellInds(i)).dc)),  '.', 'Color', colorCode(allCellInds(i),:)); hold on;
    axis ([0 2 -1.6 2.2]);
    set(gca, 'YTick', yTickPos, 'YTickLabel', [], 'XTick', []);
end
kFilt.select();
plot(-49:0, mean(kMat,2), 'k', 'LineWidth' , 2); 
hFilt.select();
plot(cellModels(1).iht, mean(ihMat,2), 'k', 'LineWidth' , 2); 
bFilt.select();
plot(1, mean(dcComp), 'k.', 'MarkerSize', 12);


currGlmPan.margin = 5;
end
glmParmPan.margin = 5;
glmParmPan.fontsize = 6;


%% compute pca space of neurons

allFilts = [kMat' downsample(ihMat(101:end,:),10)' dcComp']';
allFilts = zscore(allFilts);

data = allFilts;
allFiltsMean = mean(allFilts,2);
pcaEigs = 2;
[outputPCsA,pMatrixA, rankEigen] = princomp(data');
pctVarExpA = cumsum(rankEigen)./sum(rankEigen);

%% plot 1st 3 PCAs

% plot pct variance explained by cumulative PCAs
figure;
pctVarExpPan = panel();
pctVarExpPan.select();
plot(pctVarExpA(1:15)*100, '.')

pctVarExpPan.fontsize = 6;

pctVarExpPan.xlabel('PCA index');
pctVarExpPan.ylabel('Pct variance explained');

axis([0 16 0 100]);



figure;
pcaPlotPan = panel();
pcaPlotPan.select();
plot(outputPCsA(:,1:3))
pcaPlotPan.fontsize = 6;
pcaPlotPan.xlabel('Parameter index');
pctVarExpPan.ylabel('Parameter value (a.u.)');



% plot neurons along PCs 1&2 and 1&3
figure;
plotPCs = [1 2; 1 3];
neuronDistnPan = panel();
neuronDistnPan.pack('h', size(plotPCs, 1));

for i = 1:size(plotPCs, 1)

    currPlotPCs = plotPCs(i,:);

pMatrix = pMatrixA;
outputPCs = outputPCsA;
neuronDistnPan(i).select();

hold on;
for j = 1:size(pMatrix,1)
    plot(pMatrix(j,currPlotPCs(1)), pMatrix(j,currPlotPCs(2)), '.', 'Color', colorCode(j,:,:), 'MarkerSize', 10);
end

xlabel(['PCA ', num2str(currPlotPCs(1))]);
ylabel(['PCA ', num2str(currPlotPCs(2))]);
axis([-10 10 -4 6]);

end
neuronDistnPan.fontsize = 6;