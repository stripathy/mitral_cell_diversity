numNeurons = 40;
colorCode = colormap(jet(length(keepCells)));

figure;
plot(1:numNeurons, repmat(stimEntPerSec, 1, numNeurons),'k');

hold on;
neuronMatHomo = reshape(totalEntPerSec(1:440), 40, 11)';
for i = 1:11
    plot(1:numNeurons, neuronMatHomo(i,:), 'Color', colorCode(neuronsToSim(i),:));
end

neuronMatHetero = reshape(totalEntPerSec(441:end), 50, 40);
meanEntHetero = mean(neuronMatHetero, 1);

plot(1:numNeurons, meanEntHetero, 'r*');

numHomoNeurons = 11;
figure;
hold on;
neuronMatHomo = reshape(meanReconAcc(1:numHomoNeurons*numNeurons), numNeurons, numHomoNeurons)';
for i = 1:numHomoNeurons
    plot(1:numNeurons, neuronMatHomo(i,:), 'Color', colorCode(neuronsToSim(i),:));
end

neuronMatHetero = reshape(meanReconAcc(numHomoNeurons*numNeurons+1:end), numPopsPerSize, numNeurons);
meanEntHetero = mean(neuronMatHetero, 1);

plot(1:numNeurons, meanEntHetero, 'r*');
xlabel('population size'); ylabel('reconstruction error (rmse)');
axis([0 40 .3 .9]);


currStimInd = 26;
homoCellInd = 2;
plotInds(1) = homoCellInd*maxPopSize*numSignals - numSignals + currStimInd;
plotInds(2) = length(reconStructIn) - numSignals + currStimInd;

colorCodePlot(1,:) = colorCode(neuronsToSim(homoCellInd),:);
colorCodePlot(2,:) = [1 0 0];

figure; subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2) % axis([1 slen -3 3]);
testSpikesPlot = [];
hold on;
rasterColors = [];
for i = 1:2%length(currCellInds)
    tempCurrCellInds = reconStructIn(plotInds(i)).popMakeup;
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim,  'Color', colorCodePlot(i,:), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    
%     for j = 1:length(tempCurrCellInds)
%         rasterColors(end+1,:) = colorCode(tempCurrCellInds(j),:);
%     end
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