    
% compute info spect density between stim and decoded residuals
Fs = 1000;
nfft = 2^nextpow2(slen);
f = Fs/2*linspace(0,1,nfft/2+1);

StimFFT = fft(outSignal,nfft)/slen;
meanStimFFT = mean(2*abs(StimFFT(1:nfft/2+1,:)),2);

for i = 1:numPops
    for j = 1:numSignals
        reconStim(:,j) = reconStructOut(j+(i-1)*numSignals).optStim;
    end
    reconResids = reconStim-outSignal;
    residFFT = fft(reconResids,nfft)/slen;
    meanResidFFT = mean(2*abs(residFFT(1:nfft/2+1,:)),2);
    
    infoSpect(:,i) = log2(meanStimFFT./meanResidFFT);
end


figure;
hold on;
for i = 1:length(cellPlotInds)
    currCell = cellPlotInds(i);
    plot(histEdges, isiHist(:,currCell), 'Color', colorCode(currCell,:))
end
xlabel('time (ms)'); ylabel('ISI probability');

tapers=3; avg=1; fs = 1000;
    pars = struct ('Fs', fs, ...
        'tapers', [tapers 2*tapers-1], ...
        'pad', 0, ...
        'trialave', 1);
    
clear currStim
for j = 1:length(cellModels)
    %compute coherence for each stim and all repeats
    for k = 1:numSignals
        cellStimInd = k + (j-1)*numSignals;
        currStim(:, 1:numSimRepeats) = repmat(outSignal(:,k), 1, numSimRepeats);
        currSpikes = simSpikes(cellStimInd).dat*.001;
        spikeStruct = [];
        for i = 1:size(currSpikes,2)
            spikeStruct(i).times(:,1) = nonzeros(currSpikes(:,i));
        end
        [C(:,k),phi(:,k),S12,S1(:,k),S2(:,k),freq]=coherencycpt(currStim, spikeStruct, pars);
    end
    meanCohere(:,j) = nanmean(C');
    meanPhase(:,j) = nanmean(phi');
    meanSpect2(:,j) = mean(S2');
    meanSpect1(:,j) = mean(S1');
end


freqPts = 1:72;

figure;
hold on;
for i = 1:length(cellPlotInds)
    currCell = cellPlotInds(i);
    plot(freq(:), meanCohere(:,currCell), 'Color', colorCode(currCell,:) , 'LineWidth', 2)
end
xlabel('frequency (Hz)'); ylabel('stimulus-spiketrain coherence');

figure;
hold on;
for i = 1:length(cellPlotInds)
    currCell = cellPlotInds(i);
    plot(freq(:), meanSpect2(:,currCell), 'Color', colorCode(currCell,:) , 'LineWidth', 2)
end
xlabel('frequency (Hz)'); ylabel('spike train power');

plot(freq(:), meanSpect1(:,1), 'k', 'LineWidth', 2)
xlabel('frequency (Hz)'); ylabel('stimulus power');

figure;
hold on;
for i = 1:length(cellPlotInds)
    currCell = cellPlotInds(i);
    plot(freq(:), meanPhase(:,currCell), 'Color', colorCode(currCell,:) , 'LineWidth', 2)
end
xlabel('frequency (Hz)'); ylabel('stimulus-spiketrain phase coherence');

figure;
hold on;
for i = 1:length(cellPlotInds)
    currCell = cellPlotInds(i);
    plot(freq(:), infoSpect(:,currCell), 'Color', colorCode(currCell,:) , 'LineWidth', 2)
end
xlabel('frequency (Hz)'); ylabel('info spectra (bits/Hz)');