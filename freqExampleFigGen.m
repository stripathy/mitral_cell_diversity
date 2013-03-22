colorCodePlot(1,:) = [1 0 0];
colorCodePlot(2,:) = [0 1 0];
colorCode = colormap(jet(length(keepCells)));

stimStatExamples = [1 2 6 4 10];
numFreqs = length(stimStatExamples);
figure;
freqInds1 = 1:53;
freqInds2 = 1:257;
% freqPan = panel();
% freqPan.pack('v',3)
for k = 1:numFreqs
    currStimStat = stimStatExamples(k);        

    currFreqInd = find(stimStatExamples(k) == freqSimStimInds);
%     freqPan(k).select();
    subplot(2,numFreqs,k);
    plot(freq(freqInds1), stimSpect(currFreqInd).dat(freqInds1), 'k', 'LineWidth' , .5); hold on;
    for i = 1:2%length(currCellInds)
        plot(freq(freqInds1), reconSpect(currFreqInd,i).dat(freqInds1), 'Color', colorCodePlot(i,:), 'LineWidth' , .5);
    end
    axis tight;
    set(gca, 'FontSize', 6);
    box off;
    
    subplot(2,numFreqs,k+numFreqs);
    loglog(freq(freqInds2), stimSpect(currFreqInd).dat(freqInds2), 'k', 'LineWidth' , .5); hold on;
    for i = 1:2%length(currCellInds)
        loglog(freq(freqInds2), reconSpect(currFreqInd,i).dat(freqInds2), 'Color', colorCodePlot(i,:), 'LineWidth' , .5);
    end
    axis([0 500 .001 1]);
    set(gca, 'XTick', [1 10 100]);
    set(gca, 'YTick', [.01 .1]);
    set(gca, 'FontSize', 6);
    
    box off;
end

