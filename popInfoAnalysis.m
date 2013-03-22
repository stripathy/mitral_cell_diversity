%% compute mutual info and info spect density

% compute entropy using MAP residuals - CHECK
stimEnt = .5*log2(det(stimCovMat));
for i = 1:numPops
    for j = 1:numSignals
        inputStim = outSignal(:,j);
        reconStim = reconStructOut(j+(i-1)*numSignals).optStim;
        tempInfos(j) =  reconInfo(inputStim,  reconStim);
    end
    info(i) = mean(tempInfos);
end

a = reshape(totalEnt, 12,3);
 
% compute entropy using laplace approx on posterior MAP distn
stimEnt = .5*sum(myLog2(eig(stimCovMat))) + .5*slen*myLog2(2*pi*exp(1));
for i = 1:numPops
%     hessDetSample = zeros(1, numSignals);
    respEntSample = zeros(1, numSignals);
    for j = 1:numSignals
%         hessMatTemp = reconStructOut(j+(i-1)*numSignals).hessian;
        respEntSample(j) = -.5*sum(myLog2(reconStructOut(j+(i-1)*numSignals).hessianEigs));
        meanReconAccSample(j) = sqrt(mean((outSignal(:,j) - reconStructOut(j+(i-1)*numSignals).optStim).^2));
    end
%     [respEntCI(i,1:2)] = myBootstrapCI(respEntSample);
%     respEnt(i) = mean(-.5*myLog2(hessDetSample)) + .5*slen*myLog2(2*pi*exp(1));
    respEnt(i) = mean(respEntSample) + .5*slen*myLog2(2*pi*exp(1));
    respEntError(i) = std(respEntSample)/sqrt(numSignals);
    totalEnt(i) = stimEnt - respEnt(i);
    totalEntError(i) = respEntError(i);
    
    meanReconAcc(i) = mean(meanReconAccSample);
    meanReconAccError(i) = std(meanReconAccSample)/sqrt(numSignals);
%     totalEntCI(i,:) = abs(respEntCI(i,:) - respEnt(i));
end

stimEntPerSec = stimEnt/(slen*.001);
totalEntPerSec = totalEnt/(slen*.001);
totalEntPerSecError = totalEntError/(slen*.001);
% totalEntCI = totalEntCI/(slen*.001);


%% compute mean reliability of each cell

%% compute ISI dist and mean reliability of each cell
deltaLen = 6;
testInterval = [0 slen];
% [gammaModel(i), gammaNorm(i), gammaReal(i)] = gammaCoincFact(realSpikes, modelSpikes, testInterval, deltaLen);

isiBinSize = 2;
isiMaxTime = 100;
for j = 1:length(cellModels)
    %compute coherence for each stim and all repeats
    fullSpikeMat = [];
    for k = 1:numSignals
        cellStimInd = k + (j-1)*numSignals;
        currSpikes = simSpikes(cellStimInd).dat;
%         [blah1, blah1, tempReliab(k)] = gammaCoincFact(currSpikes, currSpikes, testInterval, deltaLen);
        fullSpikeMat = vectCat(fullSpikeMat, currSpikes);
    end
    [isiHist(:,j), histEdges] = myISI(fullSpikeMat, isiBinSize, isiMaxTime);
    meanFR(j) = nnz(fullSpikeMat)/(numSimRepeats*numSignals*slen*.001);
    
%     cellReliab(j)  = mean(tempReliab);
end

cellPlotInds = [16 32];


%compute mean zero lag correlation between neuron spikes
corrBoxSize = 5;
zeroLagCorrMatFull = zeros(numModels, numModels, numSignals);
parfor k = 1:numSignals
    %get spike trains for each neuron each trial
    spikeMat = [];
    for j = 1:length(cellModels)
        cellStimInd = k + (j-1)*numSignals;
        currSpikes = simSpikes(cellStimInd).dat;
        spikeMat = vectCat(spikeMat, currSpikes);
    end
    tempCorrMat = myCorr1(spikeMat, corrBoxSize);
    [zeroLagCorrMatFull(:,:,k)] = corrMatProcess(tempCorrMat, numModels, numSimRepeats);
end
zeroLagCorrMat = mean(zeroLagCorrMatFull,3);
corrMat = zeroLagCorrMat;
symCorrMat = corrMat + triu(corrMat,1)';

%test corrMat by plotting random spike rasters from cells
testCells = [13 14];
clear rasterColors
spikeMat = [];
k = randi(numSignals, 1, 1);
for i = 1:length(testCells)
    cellStimInd = k + (testCells(i)-1)*numSignals;
    currSpikes = simSpikes(cellStimInd).dat;
    spikeMat = vectCat(spikeMat, currSpikes);
end
figure;
createRaster2(spikeMat', 0, slen);


%compute mean pairwise gamma factor between neurons across stims
deltaLen = 6;
testInterval = [0 slen];
gammaFactMatAll = zeros(numModels, numModels, numSignals);
gammaFactMatNormAll = zeros(numModels, numModels, numSignals);

parfor k = 1:numSignals
    %get spike trains for each neuron each trial
    gammaFactMatTemp = zeros(length(cellModels), length(cellModels));
    gammaFactMatNormTemp = zeros(length(cellModels), length(cellModels));
    for j = 1:length(cellModels)
        neuron1Ind = k + (j-1)*numSignals;
        neuron1Spikes = simSpikes(neuron1Ind).dat;
        for i = 1:length(cellModels)
            neuron2Ind = k + (i-1)*numSignals;
            neuron2Spikes = simSpikes(neuron2Ind).dat;
            [gammaFactMatTemp(j,i), gammaFactMatNormTemp(j,i), blah1] = gammaCoincFact(neuron1Spikes, neuron2Spikes, testInterval, deltaLen);
        end
    end
    gammaFactMatAll(:,:,k) = gammaFactMatTemp;
    gammaFactMatNormAll(:,:,k) = gammaFactMatNormTemp;
end
gammaFactMat = mean(gammaFactMatAll,3);
gammaFactMatNorm = mean(gammaFactMatNormAll,3);


%compute mean pairwise correlation between reconstructions (how similar are
%reconstuctions for each pair of neurons?)

%not finished
plotInds = currStimInd + numSignals*(currCellInds-1);

parfor k = 1:numSignals
    %get spike trains for each neuron each trial
    reconMat = zeros(slen, numSimRepeats*numSignals);
    for j = 1:length(cellModels)
        cellStimInd = k + (j-1)*numSignals;
        tempInds = (j-1)*numSimRepeats + 1:j*numSimRepeats;
        currStim(:, tempInds) = repmat(outSignal(:,k), 1, numSimRepeats);
        currSpikes = simSpikes(cellStimInd).dat;
        spikeMat = vectCat(spikeMat, currSpikes);
    end
    tempCorrMat = myCorr1(spikeMat, corrBoxSize);
    [zeroLagCorrMatFull(:,:,k)] = corrMatProcess(tempCorrMat, numModels, numSimRepeats);
end
zeroLagCorrMat = mean(zeroLagCorrMatFull,3);

%% compute 0-lag correlation values (rho) between cells using corrFxns
corrMat = zeros(numModels, numModels);
cnt = 1;
for i = 1:numModels
    for j = 1:numModels
        if i==j
            corrMat(i,j) = computeRhos(xcorrs(i,:),acorrs1(i,:), acorrs2(i,:));
        else
            corrMat(i,j) = computeRhos(xcorrsPairwise(cnt,:) ,autoCorrs(i,:),autoCorrs(j,:));
        end
        cnt = cnt + 1;
    end
end

%% plot single cell info stats

colorCode = colormap(jet(length(keepCells)));
figure;
hold on;
for i = 1:length(colorCode)
    bar(i, totalEntPerSec(i), 'FaceColor', colorCode(i,:))
    errorbar(i, totalEntPerSec(i), totalEntPerSecCI(i,1), totalEntPerSecCI(i,2),'k');
end
xlabel('cell number'); ylabel('info (bits/sec)');


figure;
hold on;
for i = 1:length(colorCode)
    errorbar(meanFR(i), totalEntPerSec(i), totalEntPerSecCI(i,1), totalEntPerSecCI(i,2),'k');
    plot(meanFR(i), totalEntPerSec(i),'*', 'Color', colorCode(i,:), 'LineWidth', 2)
end
xlabel('mean FR (Hz)'); ylabel('info (bits/sec)');

neuronInfoPerSpike = totalEntPerSec(1:numModels)./meanFR;
neuronInfoPerSpikeCI(1:numModels,1) = totalEntPerSecCI(1:numModels,1)./meanFR';
neuronInfoPerSpikeCI(1:numModels,2) = totalEntPerSecCI(1:numModels,2)./meanFR';
% figure out a regression line for info (bits/spike) vs reliability
% [linReliabFact, bint, rVal, rint, stats] = regress(neuronInfoPerSpike', cellReliab'); 

cellReliab2 = diag(corrMat);
figure;
hold on;
for i = 1:length(colorCode)
    errorbar(cellReliab2(i), neuronInfoPerSpike(i), neuronInfoPerSpikeCI(i,1), neuronInfoPerSpikeCI(i,2),'k');
    plot(cellReliab2(i), neuronInfoPerSpike(i),'*', 'Color', colorCode(i,:), 'LineWidth', 2)
end
xlabel('neuron reliability'); ylabel('info (bits/spike)');
% plot(0:.1:.7, (0:.1:.7)*linReliabFact)


figure;
hold on;
for i = 1:length(colorCode)
    bar(i, neuronInfoPerSpike(i) - linReliabFact*cellReliab(i), 'FaceColor', colorCode(i,:))
end
xlabel('cell number'); ylabel('info (bits/sec/linear reliability)');


       
%% plot some recons and spike trains

colorInds = ['bgr'];
rasterColorInds(1) = colorInds(1);
rasterColorInds(2:3) = colorInds(2);
rasterColorInds(4:9) = colorInds(3);

myCellInd = 7;

currCellInd = find(keepCells == myCellInd);
% currCellInd = 3;
currStimInd = randi(100,1);

figure; subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2)
% axis([1 slen -3 3]);

plotInds = numSignals*numModels*[0:2] + currStimInd + numSignals*(currCellInd-1);
testSpikesPlot = [];
hold on;
for i =1:3
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim, colorInds(i), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
end
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus value (zscore)');
subplot(3,1,3);createRaster2(testSpikesPlot',1, slen, rasterColorInds); ylabel('Neuron ID');
xlabel('time (ms)');   


%% some analysis for population Sims (> 1 cell).

entMat = zeros(length(keepCells), length(keepCells));
% entMatNorm = zeros(length(keepCells), length(keepCells));
% diffEntMat = zeros(length(keepCells), length(keepCells)); %this is the differential entropy matrix (AandB - B)/A

cnt = 1;
for i = 1:length(keepCells)
    for j = i:length(keepCells)
        entMat(i,j) = totalEnt(cnt+length(keepCells));
        entMatPerSec(i,j) = totalEntPerSec(cnt+length(keepCells));
        entMatPerSec(j,i) = entMatPerSec(i,j);
        entMatPerSecError(i,j) = totalEntPerSecError(cnt);
        entMatPerSecError(j,i) = entMatPerSecError(i,j);
        synMatPerSec(i,j) = (entMatPerSec(i,j) - totalEntPerSec(j) - totalEntPerSec(i))/(entMatPerSec(i,j));
        synMatPerSec(j,i) = synMatPerSec(i,j);
        tempError = (sqrt(entMatPerSecError(i,j)^2 + totalEntPerSecError(i)^2 + totalEntPerSecError(j)^2)/(entMatPerSec(i,j) - totalEntPerSec(j) - totalEntPerSec(i)))^2 ...
            + ( entMatPerSecError(i,j)/entMatPerSec(i,j) )^2;
        synMatPerSecError(i,j) = abs(sqrt(tempError)*synMatPerSec(i,j));
        synMatPerSecError(j,i) = synMatPerSecError(i,j);
        cnt = cnt + 1;
    end
end


for i = 1:length(keepCells)
    for j = 1:length(keepCells)
%         diffEntMat(i,j) = (entMatPerSec(i,j) - totalEntPerSec(j))/totalEntPerSec(i);
        diffEntMat(i,j) = (entMatPerSec(i,j) - totalEntPerSec(j) - totalEntPerSec(i))/(totalEntPerSec(j) + totalEntPerSec(i));
%         diffEntMat(i,j) = (entMatPerSec(i,j) - totalEntPerSec(j));
    end
end
figure;
imagesc(1:numModels, 1:numModels, diffEntMat(1:numModels,1:numModels)); colorbar; 
caxis([.5 1.3]); 
title('differential information (AandB - B)/A, normalized');

figure; plot(symCorrMat(1:end, 1:end), diffEntMat(1:end, 1:end),'.')
figure; plot(gammaFactMat(3:end, 3:end), diffEntMat(3:end, 3:end),'.')


figure;
imagesc(1:numModels, 1:numModels, entMatPerSec); colorbar; title('pairwise mutual info (bits/sec)');

figure;
imagesc(1:32, 1:32, symCorrMat(1:end, 1:end)); colorbar; title('spike time correlation matrix (rho)');

figure;
imagesc(entMatNorm); colorbar; title('pairwise mutual info (bits/spike)');

figure;
imagesc(entMatNormNorm); colorbar; title('pairwise mutual info norm (bits/spike)');

figure;
imagesc(3:32, 3:32, entMatNormReliab(3:end, 3:end)); colorbar; title('pairwise mutual info (bits/spike/per unit reliability)');

figure;
imagesc(entMatPerSecSum-entMatPerSec); colorbar; title('pairwise mutual info A+B - AandB (bits/sec)');

figure;
imagesc(entMatPerNormSum-entMatNorm); colorbar; title('pairwise mutual info A+B - AandB (bits/spike)');

figure;
imagesc((entMatPerNormSum-entMatNorm)./entMatPerNormSum); colorbar; title('pairwise mutual info A+B - AandB (bits/spike, normalized)');

figure;
plot(nonzeros(corrMat), nonzeros(entMatNormNorm),'.'); title('normalized pairwise information gain vs correlation');

figure;
plot(nonzeros(corrMat), nonzeros(entMatPerNormSum-entMatNorm),'.')

entMatPerSec = entMatPerSec + triu(entMatPerSec,1)';




%% plot some final figs
figure;
hold on;
for i = 1:length(colorCode)
    plot(cellReliab(i), diffEntMat(i,i),'*', 'Color', colorCode(i,:), 'LineWidth', 2)
end
xlabel('neuron reliability'); ylabel('differential information'); %axis([0 1 0 1.5]);

%% plotting spikes and recons for populations

% plot a histogram of diffEntVals for a given cell
currCell = 35;
compCell = 35;

% figure; hist(diffEntMat(currCell, :),15); xlabel('differential information'); axis tight
% figure; plot(symCorrMat(currCell, :),diffEntMat(currCell, :),'.'); 
% xlabel('correlation'); ylabel('differential information');


%plot pairs of rasters from diff cells
%test corrMat by plotting random spike rasters from cells

corrTpts = -50:50;
% figure; plot(corrTpts, xcorrs(currCell,:)', corrTpts, xcorrsPairwise((currCell-1)*numModels + compCell,:)');
figure; plot(corrTpts, xcorrs(currCell,:)', corrTpts, fliplr(xcorrsPairwise((compCell-1)*numModels + currCell,:))');


testCells = [currCell compCell];
clear rasterColors
rasterColors(1:numSimRepeats) = 'b';
rasterColors(numSimRepeats+1:2*numSimRepeats) = 'g';
spikeMat = [];
currStimInd = randi(numSignals, 1, 1);
for i = 1:length(testCells)
    cellStimInd = currStimInd + (testCells(i)-1)*numSignals;
    currSpikes = simSpikes(cellStimInd).dat;
    spikeMat = vectCat(spikeMat, currSpikes);
end
figure;
createRaster2(spikeMat', 0, slen, rasterColors);


popIndMat = zeros(numModels, numModels);
cnt = 1;
for i = 1:numModels
    for j = i:numModels
        popIndMat(i,j) = cnt + numModels;
        popIndMat(j,i) = cnt + numModels;
        cnt = cnt + 1;
    end
end
colorCodePlot = [0 0 1; 1 0 0; 0 1 0];

% currStimInd = randi(numSignals, 1, 1);
% currCellInds = [currCell 0; currCell testCells(2)];
currCellInds = [currCell 0; compCell 0; currCell testCells(2)];
% currCellInd = 3;
 

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
for i = 1:length(currCellInds)
    tempCurrCellInds = currCellInds(i,:);
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim, 'Color', colorCodePlot(i,:), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    testSpikesPlot = vectCat(testSpikesPlot, 0);
    rasterColors(end+1:end+nnz(tempCurrCellInds),:) = repmat(colorCodePlot(i,:), nnz(tempCurrCellInds),1);
    rasterColors(end+1,:) = [1,1,1]; %white
%     rasterColors(end+1:end+nnz(tempCurrCellInds),:) = colorCode(nonzeros(tempCurrCellInds),:);
%     rasterColors(end+1,:) = [1,1,1]; %white
end
axis([0 slen -3 3]);
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus value (zscore)');
subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('time (ms)');



% plot GLM fit parms
figure; 
allCellInds = unique(nonzeros(currCellInds));
for i = 1:length(allCellInds)
    
    subplot(1,3,1);plot(-49:0, cellModels(allCellInds(i)).k, 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 2); axis([-50 0 -.5 1.5]); hold on;
    subplot(1,3,2);plot(cellModels(1).iht,cellModels(allCellInds(i)).ihbas*cellModels(allCellInds(i)).ih, 'Color', colorCode(allCellInds(i),:), 'LineWidth' , 2); axis([0 63 -5 3.2]); hold on;
    subplot(1,3,3);plot(1, cellModels(allCellInds(i)).dc,  '*', 'Color', colorCode(allCellInds(i),:)); hold on;
end

%% generate all pairwise crosscorrelations

% This code assumes that there is a single long stimulus and multiple (even
% number) of repeats

%get all spikes from all neurons
allSpikes = [];
for i = 1:length(keepCells)
    allSpikes = vectCat(allSpikes, simSpikes(i).dat);
end
trialsPerNeuron(1:numModels) = numSimRepeats;

% run xcorr algorithm
tic
[xcorrs, acorrs1, acorrs2, xcorrsPairwise, autoCorrs]=crossCorrAllCells(allSpikes, trialsPerNeuron, slen);
toc



%% play with correlogram functions
tempSpikes = simSpikes(9).dat;  
tempSpikes2 = simSpikes(9).dat(:, randi(numSimRepeats,numSimRepeats,1));

% boxSize = 8;
% [corr] = myCorr1(tempSpikes, boxSize);

tmax = slen;
[ac1,ac2,xc,ac1all,ac2all,xcall]=crossCorrSonya(tempSpikes',tempSpikes2',tmax);
plot(ac1)

tmax = slen;
[ac1,ac2,xc,ac1all,ac2all,xcall]=crossCorrSonyaTest(tempSpikes',tempSpikes2',tmax);
plot(ac1)



tempStim = outSignal(:,1);
staWin = 50;
[filtOut, offset] = myCorrelogram(tempSpikes, endPts, staWin)

spikeList = sort(nonzeros(tempSpikes));

tempSpikesCorr = zeros(2, 16, 50);
for i = 1:50
    tempSpikesCorr(1,1:length(tempSpikes(:,i)),i) = tempSpikes(:,i)';
    tempSpikesCorr(2,1:length(tempSpikes2(:,i)),i) = tempSpikes2(:,i)';
end

tempSpikes = simSpikes(2205).dat;
a = myCorr1(tempSpikes);
imagesc(a)
mean(nonzeros(triu(a)))

[rhos,acorrs,xcorrs,Instrate] = correlateTest(tempSpikesCorr,256, [50 50]);

[filtOut, offset] = myCorrelogram(spikes, endPts, staWin);

%% compute total entropy but select for only sims with 
% compute entropy using laplace approx on posterior MAP distn
stimEnt = .5*sum(myLog2(eig(stimCovMat))) + .5*slen*myLog2(2*pi*exp(1));
for i = 1:numPops
%     hessDetSample = zeros(1, numSignals);
    respEntSample = zeros(1, numSignals);
    for j = 1:numSignals
%         hessMatTemp = reconStructOut(j+(i-1)*numSignals).hessian;
        % check if pop has spikes or not, if not, continue
        if sum(sum(reconStructIn(j+(i-1)*numSignals).testSpikes) == 0) == 0
            respEntSample(j) = -.5*sum(myLog2(reconStructOut(j+(i-1)*numSignals).hessianEigs));
            meanReconAccSample(j) = sqrt(mean((outSignal(:,j) - reconStructOut(j+(i-1)*numSignals).optStim).^2));
            %c = 'ab'
        else
            i, j
            respEntSample(j) = NaN;
            meanReconAccSample(j) = NaN;
        end
    end
%     [respEntCI(i,1:2)] = myBootstrapCI(respEntSample);
%     respEnt(i) = mean(-.5*myLog2(hessDetSample)) + .5*slen*myLog2(2*pi*exp(1));
    respEnt(i) = nanmean(respEntSample) + .5*slen*myLog2(2*pi*exp(1));
    respEntError(i) = nanstd(respEntSample)/sqrt(numSignals);
    totalEnt(i) = stimEnt - respEnt(i);
    totalEntError(i) = respEntError(i);
    
    meanReconAcc(i) = nanmean(meanReconAccSample);
    meanReconAccError(i) = nanstd(meanReconAccSample)/sqrt(numSignals);
%     totalEntCI(i,:) = abs(respEntCI(i,:) - respEnt(i));
end

stimEntPerSec = stimEnt/(slen*.001);
totalEntPerSec = totalEnt/(slen*.001);
totalEntPerSecError = totalEntError/(slen*.001);

%% stuff to compute CV ISI

for i = 1:length(cellModels)
    spikeMat = simSpikes(i).dat;
    isis = diff(spikeMat, 1, 1);
    isis = isis(isis>0);
    cvisi(i) = std(isis)/mean(isis);
    
end

colorCode = colormap(jet(length(keepCells)));
figure;
hold on;
% errorbar(meanFR, totalEntPerSec(1:numModels), totalEntPerSecError(1:numModels),'.k');
for i = 1:numModels
    plot(totalEntPerSec(i), cvisi(i), '*', 'Color', colorCode(i,:), 'LineWidth', 2)
end
xlabel('info (bits/sec)'); ylabel('CV ISI distribution'); 

colorCode = colormap(jet(length(keepCells)));
figure;
hold on;
% errorbar(meanFR, totalEntPerSec(1:numModels), totalEntPerSecError(1:numModels),'.k');
for i = 1:numModels
    plot(cvisi(i), totalEntPerSec(i), '*', 'Color', colorCode(i,:), 'LineWidth', 2)
end
ylabel('info (bits/sec)'); xlabel('CV ISI distribution'); 