% plot figures based on population sims
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\pairwiseInfo2\allPairwisePopsFull.mat'
%% plot figures for single cell infos

% info by FR
colorCode = colormap(jet(length(keepCells)));
figure;
hold on;
errorbar(meanFR, totalEntPerSec(1:numModels), totalEntPerSecError(1:numModels),'.k');
for i = 1:numModels
    plot(meanFR(i), totalEntPerSec(i),'*', 'Color', colorCode(i,:), 'LineWidth', 2)
end
xlabel('mean FR (Hz)'); ylabel('info (bits/sec)'); axis([0 70 0 100]);

neuronInfoPerSpike = totalEntPerSec(1:numModels)./meanFR;
neuronInfoPerSpikeError = totalEntPerSecError(1:numModels)./meanFR;
% neuronInfoPerSpikeCI(1:numModels,1) = totalEntPerSecCI(1:numModels,1)./meanFR';
% neuronInfoPerSpikeCI(1:numModels,2) = totalEntPerSecCI(1:numModels,2)./meanFR';

% info (bits/spike) by reliability
cellReliab2 = diag(corrMat);
figure;
hold on;
errorbar(cellReliab2, neuronInfoPerSpike, neuronInfoPerSpikeError,'.k');
for i = 1:numModels
    plot(cellReliab2(i), neuronInfoPerSpike(i),'*', 'Color', colorCode(i,:), 'LineWidth', 2)
end
xlabel('neuron reliability'); ylabel('info (bits/spike)'); axis([0 .5 .5 3]);
% plot(0:.1:.7, (0:.1:.7)*linReliabFact)

% plot a stimulus and 2 example recons for two cells

currCell = 18;
compCell = 35;

% figure; hist(diffEntMat(currCell, :),15); xlabel('differential information'); axis tight
% figure; plot(symCorrMat(currCell, :),diffEntMat(currCell, :),'.'); 
% xlabel('correlation'); ylabel('differential information');


%plot pairs of rasters from diff cells
%test corrMat by plotting random spike rasters from cells

popIndMat = zeros(numModels, numModels);
cnt = 1;
for i = 1:numModels
    for j = i:numModels
        popIndMat(i,j) = cnt + numModels;
        popIndMat(j,i) = cnt + numModels;
        cnt = cnt + 1;
    end
end
colorCodePlot(1,:) = colorCode(currCell,:);
colorCodePlot(2,:) = colorCode(compCell,:);

corrTpts = -50:50;
% figure; plot(corrTpts, xcorrs(currCell,:)', corrTpts, xcorrsPairwise((currCell-1)*numModels + compCell,:)');
figure; hold on;
for i = 1:2
    if i == 1
        tempCell = currCell;
    else
        tempCell = compCell;
    end
 plot(corrTpts, xcorrs(tempCell,:), 'Color', colorCodePlot(i,:), 'LineWidth', 2)
end
xlabel('time (ms)'); ylabel('spike probability'); axis tight;

currStimInd = 12;%randi(numSignals, 1, 1);


currCellInds = [currCell 0; compCell 0];

for i = 1:size(currCellInds,1)
    tempCurrCellInds = currCellInds(i,:);
    if nnz(tempCurrCellInds) < 2
        plotInds(i) = currStimInd + numSignals*(tempCurrCellInds(1)-1);
    else %two cells
        plotInds(i) = currStimInd + numSignals*(popIndMat(tempCurrCellInds(1), tempCurrCellInds(2))-1);
    end
end 

figure; subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2) % axis([1 slen -3 3]);
testSpikesPlot = [];
hold on;
rasterColors = [];
for i = 1:2%length(currCellInds)
    tempCurrCellInds = currCellInds(i,:);
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim, 'Color', colorCodePlot(i,:), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    
    rasterColors(end+1:end+nnz(tempCurrCellInds),:) = repmat(colorCodePlot(i,:), nnz(tempCurrCellInds),1);
    if i ~= 2
        rasterColors(end+1,:) = [1,1,1]; %white
        testSpikesPlot = vectCat(testSpikesPlot, 0);
    end
end
axis([0 slen -3 3]);
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus (a.u.)'); set(gca,'XTickLabel', []); axis tight; box off;
subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('time (ms)'); set(gca,'YTick', []);

figure; 
allCellInds = unique(nonzeros(currCellInds));
for i = 1:length(allCellInds)
    
    subplot(1,5,1:2);plot(-49:0, cellModels(allCellInds(i)).k, 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 2); axis([-50 0 -.5 1]); hold on;
    subplot(1,5,3:4);plot(cellModels(1).iht,log10(exp(cellModels(allCellInds(i)).ihbas*cellModels(allCellInds(i)).ih)), 'Color', colorCode(allCellInds(i),:), 'LineWidth' , 2); hold on;
    axis ([0 60 -1.6 2.2]);
    yTickPos = [-1:2]; yTickLabels = [.1 1 10 100];
    set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabels);
    subplot(1,5,5);plot(1, log10(exp(cellModels(allCellInds(i)).dc)),  '*', 'Color', colorCode(allCellInds(i),:)); hold on;
    axis ([0 2 -1.6 2.2]);
    set(gca, 'YTick', yTickPos, 'YTickLabel', [], 'XTick', []);
end

%% single cell infos, how does k/h influence info?
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\stimHistRatio.mat';

% figure out new color code based on sorted k/h values
colorCode = colormap(jet(length(keepCells)));

[blah sortedStimHistInds] = sort(stimHistRatio);
colorCodeStimHist = colorCode(sortedStimHistInds,:);

% % color cells directly by FR
colorCodeTotal = colormap(jet(max(floor(stimHistRatio*10)) - min(floor(stimHistRatio*10)) + 1));
for i = 1:length(keepCells)
    colorCodeStimHist(i,:) = colorCodeTotal(floor(stimHistRatio(i)*10) - min(floor(stimHistRatio*10)) +1,:);
end

% info by FR
figure;
hold on;
errorbar(meanFR, totalEntPerSec(1:numModels), totalEntPerSecError(1:numModels),'.k');
for i = 1:numModels
    plot(meanFR(i), totalEntPerSec(i),'*', 'Color', colorCodeStimHist(i,:), 'LineWidth', 2)
end
xlabel('mean FR (Hz)'); ylabel('info (bits/sec)'); axis([0 70 0 100]);

neuronInfoPerSpike = totalEntPerSec(1:numModels)./meanFR;
neuronInfoPerSpikeError = totalEntPerSecError(1:numModels)./meanFR;
% neuronInfoPerSpikeCI(1:numModels,1) = totalEntPerSecCI(1:numModels,1)./meanFR';
% neuronInfoPerSpikeCI(1:numModels,2) = totalEntPerSecCI(1:numModels,2)./meanFR';

% what is the Rsquared value between bits/sec and FR?
[b,bint,r,rint,stats] = regress(totalEntPerSec(1:numModels)', [ones(44,1) meanFR']);
stats(1)

% info (bits/spike) by reliability
cellReliab2 = diag(corrMat);
figure;
hold on;
errorbar(cellReliab2, neuronInfoPerSpike, neuronInfoPerSpikeError,'.k');
for i = 1:numModels
    plot(cellReliab2(i), neuronInfoPerSpike(i),'*', 'Color', colorCodeStimHist(i,:), 'LineWidth', 2)
end
xlabel('neuron reliability'); ylabel('info (bits/spike)');
colorbar; caxis([0 2.4]);
colorLabel = 'intrinsic input / stimulus input';
h= colorbar;
set(get(h,'ylabel'),'String', colorLabel, 'Rotation', 270, 'VerticalAlignment', 'Bottom');


%% plot figures for 2 copies of the same cell


figure;
hold on;

for i = 1:numModels
    bar(i, entMatPerSec(i,i), 'FaceColor', colorCode(i,:))
    bar(i, totalEntPerSec(i), 'FaceColor', colorCode(i,:), 'EdgeColor', 'k', 'LineWidth',1.5)
%     errorbar(i, entMatPerSec(i,i), entMatPerSecError(i,i),'k');
end
errorbar(keepCells, diag(entMatPerSec), diag(entMatPerSecError), 'Color', 'k', 'LineStyle', 'none');
errorbar(keepCells, totalEntPerSec(keepCells), totalEntPerSecError(keepCells), 'Color', 'k', 'LineStyle', 'none');
xlabel('cell number'); ylabel('info (bits/sec)'); axis tight;

% how much more info for 2 cells vs 1 cell?
a = diag(entMatPerSec)./totalEntPerSec(1:numModels)';

cellReliab2 = diag(corrMat);
% [linReliabFact, bint, rVal, rint, stats] = regress(diag(synMatPerSec), [ones(44,1) cellReliab2]); 

% write code to generate regression line based on weighting errors
xMat(keepCells, 1) = 1;
xMat(keepCells, 2) = cellReliab2;
diagErrors = diag(synMatPerSecError);
for i = keepCells
    wMat(i,i) = 1/(diagErrors(i)^2);
end
bCoeff = (xMat'*wMat*xMat)\xMat'*wMat*diag(synMatPerSec);
% bCoeff = (xMat'*xMat)\xMat'*diag(synMatPerSec);
yNew = xMat*bCoeff;

figure;
hold on;
errorbar(cellReliab2, diag(synMatPerSec), diag(synMatPerSecError), 'Color', 'k', 'LineStyle', 'none');
for i = 1:numModels
%     errorbar(cellReliab2(i), synMatPerSec(i,i), synMatPerSecError(i,i),'k');
    plot(cellReliab2(i), synMatPerSec(i,i),'*', 'Color', colorCode(i,:), 'LineWidth', 2)
end

xlabel('neuron reliability'); ylabel('synergy/redundancy');
plot(0:.1:.5, bCoeff(1) + (0:.1:.5)*bCoeff(2), 'k') 
% plot(0:.1:.5, linReliabFact(1) + (0:.1:.5)*linReliabFact(2), 'k')
% 

figure;
plot(meanFR, cellReliab2,'.')


% plot a stimulus and 2 example recons: 1 single cell, 2 copies same cell

currCell = 18;
compCell = 18;

% figure; hist(diffEntMat(currCell, :),15); xlabel('differential information'); axis tight
% figure; plot(symCorrMat(currCell, :),diffEntMat(currCell, :),'.'); 
% xlabel('correlation'); ylabel('differential information');


%plot pairs of rasters from diff cells
%test corrMat by plotting random spike rasters from cells

popIndMat = zeros(numModels, numModels);
cnt = 1;
for i = 1:numModels
    for j = i:numModels
        popIndMat(i,j) = cnt + numModels;
        popIndMat(j,i) = cnt + numModels;
        cnt = cnt + 1;
    end
end
colorCodePlot(1,:) = colorCode(currCell,:);
colorCodePlot(2,:) = colorCode(compCell,:);
symbol{1} = '-';
symbol{2} = '--';


% currStimInd = randi(numSignals, 1,1);
currStimInd = 12;

currCellInds = [currCell 0; currCell compCell];

for i = 1:size(currCellInds,1)
    tempCurrCellInds = currCellInds(i,:);
    if nnz(tempCurrCellInds) < 2
        plotInds(i) = currStimInd + numSignals*(tempCurrCellInds(1)-1);
    else %two cells
        plotInds(i) = currStimInd + numSignals*(popIndMat(tempCurrCellInds(1), tempCurrCellInds(2))-1);
    end
end 

figure; subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2) % axis([1 slen -3 3]);
testSpikesPlot = [];
hold on;
rasterColors = [];
for i = 1:2%length(currCellInds)
    tempCurrCellInds = currCellInds(i,:);
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim, symbol{i}, 'Color', colorCodePlot(i,:), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    
    rasterColors(end+1:end+nnz(tempCurrCellInds),:) = repmat(colorCodePlot(i,:), nnz(tempCurrCellInds),1);
    if i ~= 2
        rasterColors(end+1,:) = [1,1,1]; %white
        testSpikesPlot = vectCat(testSpikesPlot, 0);
    end
end
axis([0 slen -3 3]);
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus (a.u.)'); set(gca,'XTickLabel', []); axis tight; box off;
subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('time (ms)'); set(gca,'YTick', []);

% plot recons with error bars for neurons 18 and 35 and stims 12 and 1
% respectively

figure; subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2) % axis([1 slen -3 3]);
testSpikesPlot = [];
hold on;
rasterColors = [];

    
for i = 1:2%length(currCellInds)
    if currCell == 18 && i == 1
        colorInd(1,:) = colorCodePlot(1,:);
        eBarInd = 1;
    elseif currCell == 18 && i == 2
        colorInd(2,:) = [0 0 1];
        eBarInd = 2;
    elseif currCell == 35 && i == 1
        eBarInd = 3;
        colorInd(1,:) = colorCodePlot(1,:);
    elseif currCell == 35 && i == 2
        eBarInd = 4;
        colorInd(2,:) = [1 0 0];
    end
    tempCurrCellInds = currCellInds(i,:);
    
    topBox = reconStructOut(plotInds(i)).optStim + reconStructOutEBars(eBarInd).stimErrorBars;
    bottomBox = reconStructOut(plotInds(i)).optStim - reconStructOutEBars(eBarInd).stimErrorBars;
    
    plot(1:slen, topBox, 1:slen, bottomBox, [1 1], [topBox(1) bottomBox(1)], [slen slen], [topBox(end) bottomBox(end)], 'Color', colorInd(i,:), 'LineWidth' , .5);
    plot(1:slen, reconStructOut(plotInds(i)).optStim, 'Color', colorInd(i,:), 'LineWidth' , 2); 
    
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    rasterColors(end+1:end+nnz(tempCurrCellInds),:) = repmat(colorInd(i,:), nnz(tempCurrCellInds),1);
    if i ~= 2
        rasterColors(end+1,:) = [1,1,1]; %white
        testSpikesPlot = vectCat(testSpikesPlot, 0);
    end
end
axis([0 slen -3 3]);
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus (a.u.)'); set(gca,'XTickLabel', []); axis tight; box off;
subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('time (ms)'); set(gca,'YTick', []);
    
%% plot figures for all pairwise infos between cells

figure;
imagesc(1:numModels, 1:numModels, entMatPerSec); 
xlabel('Cell number'); ylabel('Cell number'); axis tight;
colorLabel = 'Info (bits/sec)';
h= colorbar;
set(get(h,'ylabel'),'String', colorLabel, 'Rotation', 270, 'VerticalAlignment', 'Bottom');

figure;
imagesc(1:numModels, 1:numModels, -synMatPerSec);
xlabel('Cell number'); ylabel('Cell number'); axis tight;
colorLabel = 'Redundancy';
h= colorbar;
set(get(h,'ylabel'),'String', colorLabel, 'Rotation', 270, 'VerticalAlignment', 'Bottom');

% ask how much less redundant hetero is relative to homo
% question 1: for which cells is hetero < homo (significantly less)
pThresh = .01;

synMatPerSecError = abs(synMatPerSecError);
probMat1 = zeros(numModels, numModels);
numRands = 100;
for i = 1:numModels
    homoMean = synMatPerSec(i,i);
    homoStd = synMatPerSecError(i,i);
%     homoRands = normrnd(homoMean, homoStd, [numRands 1]);
    for j = 1:numModels
        heteroMean = synMatPerSec(i,j);
        heteroStd = synMatPerSecError(i,j);
%         heteroRands = normrnd(heteroMean, heteroStd, [numRands 1]);
%         [hyp, pVal] = ttest(heteroRands, homoRands, .01, 'left');
        erfcInput = (homoMean - heteroMean)/(sqrt(homoStd^2 + heteroStd^2));
        pVal = 1 - erfc(erfcInput)/2;
        probMat1(i,j) = pVal;
        if pVal > 1 - pThresh
            hyp = 1;
        else
            hyp = 0;
        end
        hMat1(i,j) = hyp;
    end
end
homoMinCells = find(sum(hMat1,2)<=2);
heteroMinCells = setdiff(keepCells, homoMinCells);


% question 2: for which cells is homo < hetero (significantly more)
probMat2 = zeros(numModels, numModels);
numRands = 10000;
for i = 1:numModels
    homoMean = synMatPerSec(i,i);
    homoStd = synMatPerSecError(i,i);
    homoRands = normrnd(homoMean, homoStd, [numRands 1]);
    for j = 1:numModels
        heteroMean = synMatPerSec(i,j);
        heteroStd = synMatPerSecError(i,j);
        heteroRands = normrnd(heteroMean, heteroStd, [numRands 1]);
%         [hyp, pVal] = ttest(heteroRands, homoRands, .01, 'right');
        erfcInput = (heteroMean - homoMean)/(sqrt(homoStd^2 + heteroStd^2));
        pVal = 1 - erfc(erfcInput)/2;
        probMat2(i,j) = pVal;
    end
end

figure; imagesc(probMat2); colorbar; 
xlabel('cell number'); ylabel('cell number'); axis tight;
colorLabel = 'p value';
h= colorbar;
set(get(h,'ylabel'),'String', colorLabel, 'Rotation', 270, 'VerticalAlignment', 'Bottom');
% perform tests on question 2 probMat2

% for which neurons is homo < hetero
cnt = 0;
for i = 1:numModels
    if(sum(probMat2(i,:) < .5) < 1)
        cnt = cnt + 1;
        homoWorstNeurons(cnt) = i;
    end
end
     
% for which neurons is homo < (significantly) hetero
cnt = 0;
for i = 1:numModels
    if(sum(probMat2(i,:) < .95 ) < 2)
        cnt = cnt + 1;
        homoSigWorstNeurons(cnt) = i;
    end
end

% for which neurons is hetero < (significantly) homo
cnt = 0;
pValThresh = .001;
for i = 1:numModels
    if(sum(probMat2(i,:) < pValThresh) < 1)
        cnt = cnt + 1;
        homoWorstNeurons1(cnt) = i;
    end
end


% whats the mean redundancy of homo and hetero?
homoSynElements = diag(synMatPerSec);
homoMeanSyn = mean(diag(synMatPerSec));
homoStdErrorSyn = std(diag(synMatPerSec))/sqrt(length(diag(synMatPerSec)));
heteroSynElements = [nonzeros(triu(synMatPerSec,1));  nonzeros(tril(synMatPerSec,-1))];
heteroMeanSyn = mean(heteroSynElements);
heteroStdErrorSyn = std(heteroSynElements)/sqrt(length(heteroSynElements));

figure;
errorbar(-[homoMeanSyn heteroMeanSyn], [homoStdErrorSyn heteroStdErrorSyn], 'Color', 'k', 'LineStyle', 'none');
hold on;
bar(-[homoMeanSyn heteroMeanSyn]);
ylabel('Redundancy');
currxTickLabel = {'Homog-', 'Heterog-'};
set(gca,'XTick', [1 2], 'XTickLabel', currxTickLabel)
[pval, h] = ranksum(homoSynElements, heteroSynElements)


[heteroProbDistn, heteroXes] = hist(heteroSynElements);%/numel(heteroSynElements);
heteroProbDistn = heteroProbDistn/numel(heteroSynElements);

[homoProbDistn, homoXes] = hist(diag(synMatPerSec));
homoProbDistn = homoProbDistn/numel(diag(synMatPerSec));
figure;
hold on;
bar([homoXes heteroXes], [homoProbDistn heteroProbDistn])

numHistBins = 20;
clear binEdges binCenters counts
% for i = 1:numModels
    i = 1;
    binEdges(i,:) = linspace(min(diag(synMatPerSec)),max(diag(synMatPerSec))+.001, numHistBins);
    binCenters(i,:) = (binEdges(i,1:end-1) + binEdges(i,2:end))/2;
    [counts(i,:)] = ndHistc(diag(synMatPerSec),binEdges(i,:))/numel(diag(synMatPerSec));
% end
    i = 2;
    binEdges(i,:) = linspace(min(heteroSynElements),max(heteroSynElements)+.001, numHistBins);
    binCenters(i,:) = (binEdges(i,1:end-1) + binEdges(i,2:end))/2;
    [counts(i,:)] = ndHistc(heteroSynElements,binEdges(i,:))/numel(heteroSynElements);

orderColor = ['b', 'g', 'r'];
figure; hold on;
for i = 2:-1:1
    plot(binCenters(i,:), counts(i,:), orderColor(i));
end
xlabel('Redundancy'); ylabel('counts');
title('train distances');

for i = 1:numModels
    heteroSynMatRowMean(i) = mean( setxor(synMatPerSec(i,:), synMatPerSec(i,i)));
end

figure;
hold on;
for i = 1:numModels
    plot([1 2], [synMatPerSec(i,i); heteroSynMatRowMean(i)],'*', 'Color', colorCode(i,:), 'LineWidth', 2);
    plot([1 2], [synMatPerSec(i,i); heteroSynMatRowMean(i)], '-k');
%     plot(corrMat(currCell, i),synMatPerSec(currCell, i),'*', 'Color', colorCode(i,:), 'LineWidth', 2);
end
axis([.9 2.1 -.3 .15]);
ylabel('Redundancy'); currxTickLabel = {'Homogeneous', 'Heterogeneous'};
set(gca,'XTick', [1 2], 'XTickLabel', currxTickLabel)
% 
currCell = 18;
compCell = 18;

figure;
hold on;
errorbar(corrMat(currCell, :), synMatPerSec(currCell, :), synMatPerSecError(currCell, :), 'Color', 'k', 'LineStyle', 'none');
for i = 1:numModels
    plot(corrMat(currCell, i),synMatPerSec(currCell, i),'*', 'Color', colorCode(i,:), 'LineWidth', 2);
end

xlabel('correlation'); ylabel('redundancy/synergy'); axis([-.05 .4 -.3 0]);

currCell = 35;
compCell = 35;

figure;
hold on;
for i = 1:numModels
    plot(corrMat(currCell, i),synMatPerSec(currCell, i),'*', 'Color', colorCode(i,:), 'LineWidth', 2);
end
xlabel('correlation'); ylabel('synergy/redundancy'); axis([0 max(corrMat(currCell,:))+.02 -.25 0]);


currCell = 18;
compCell = 22;

% figure; hist(diffEntMat(currCell, :),15); xlabel('differential information'); axis tight
% figure; plot(symCorrMat(currCell, :),diffEntMat(currCell, :),'.'); 
% xlabel('correlation'); ylabel('differential information');


%plot pairs of rasters from diff cells
%test corrMat by plotting random spike rasters from cells


colorCodePlot(1,:) = colorCode(currCell,:);
colorCodePlot(2,:) = colorCode(compCell,:);

figure;
hold on;
for i = 1:2
    if i == 1
        tempCell = currCell;
        plot(corrTpts, xcorrs(currCell,:)', 'Color', colorCodePlot(i,:), 'LineWidth', 2);
    else
        tempCell = compCell;
        plot(corrTpts, xcorrsPairwise((currCell-1)*numModels + compCell,:)', 'Color', colorCodePlot(i,:), 'LineWidth', 2);
    end
end
axis tight;

currStimInd = 12;%randi(numSignals, 1, 1);


currCellInds = [currCell currCell; currCell compCell];

for i = 1:size(currCellInds,1) 
    tempCurrCellInds = currCellInds(i,:);
    if nnz(tempCurrCellInds) < 2
        plotInds(i) = currStimInd + numSignals*(tempCurrCellInds(1)-1);
    else %two cells
        plotInds(i) = currStimInd + numSignals*(popIndMat(tempCurrCellInds(1), tempCurrCellInds(2))-1);
    end
end 

figure; subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2) % axis([1 slen -3 3]);
testSpikesPlot = [];
hold on;
rasterColors = [];
for i = 1:2%length(currCellInds)
    tempCurrCellInds = currCellInds(i,:);
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim, symbol{i}, 'Color', colorCodePlot(i,:), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    
    rasterColors(end+1:end+nnz(tempCurrCellInds),:) = repmat(colorCodePlot(i,:), nnz(tempCurrCellInds),1);
    if i ~= 2
        rasterColors(end+1,:) = [1,1,1]; %white
        testSpikesPlot = vectCat(testSpikesPlot, 0);
    end
end
axis([0 slen -3 3]);
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus (a.u.)'); set(gca,'XTickLabel', []); axis tight; box off;
subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('time (ms)'); set(gca,'YTick', []);

figure; 
allCellInds = unique(nonzeros(currCellInds));
for i = 1:length(allCellInds)
    
    subplot(1,5,1:2);plot(-49:0, cellModels(allCellInds(i)).k, 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 2); axis([-50 0 -.5 1]); hold on;
    subplot(1,5,3:4);plot(cellModels(1).iht,log10(exp(cellModels(allCellInds(i)).ihbas*cellModels(allCellInds(i)).ih)), 'Color', colorCode(allCellInds(i),:), 'LineWidth' , 2); hold on;
    axis ([0 60 -1.6 2.2]);
    yTickPos = [-1:2]; yTickLabels = [.1 1 10 100];
    set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabels);
    subplot(1,5,5);plot(1, log10(exp(cellModels(allCellInds(i)).dc)),  '*', 'Color', colorCode(allCellInds(i),:)); hold on;
    axis ([0 2 -1.6 2.2]);
    set(gca, 'YTick', yTickPos, 'YTickLabel', [], 'XTick', []);
end
