
fs = 1000;

tapers=5; avg=1;
    pars = struct ('Fs', fs, ...
        'tapers', [tapers 2*tapers-1], ...
        'pad', 0, ...
        'trialave', 1);
    
    segsize=1;
    [s,f]= mtspectrumc(testStim, pars);
    plot(f, s)
    
        [s,f]=psd(testStim,2048,fs,bartlett(2048),1024);
    figure; plot(f(1:200),s(1:200)); title ('original psd');
    
%     testSpikes = .001*nonzeros(reconStructIn(1).testSpikes);
testStim = outSignal(:,1);
testStim = repmat(testStim, 1, 10);
testSpikes =simSpikes(cellStimInd).dat;
testSpikes = .001*testSpikes;
segSize = .1;
for i = 1:size(testSpikes,2)
    spikeStruct(i).times(:,1) = nonzeros(testSpikes(:,i))';
end
[C,phi,S12,S1,S2,f, zerosp]=coherencycpt(testStim, spikeStruct, segSize, pars);



[C,phi,S12,S1,S2,f]=coherencycpt(testStim(:,1), nonzeros(testSpikes(:,1)), pars);
    

clear currStim
for j = 1:length(cellModels)
    %compute coherence for each stim and all repeats
    for k = 1:numSignals
        cellStimInd = k + (j-1)*numSignals;
        currStim(:, 1:numSimRepeats) = repmat(outSignal(:,k), 1, numSimRepeats);
        currSpikes = simSpikes(cellStimInd).dat*.001;
        clear spikeStruct
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
plot(freq(freqPts), meanCohere(freqPts,[10 35]))
plot(freq(freqPts), meanCohere(freqPts,:)); xlabel('frequency (Hz)'); ylabel('coherence');
plot(freq(freqPts), meanPhase(freqPts,[4 5])); xlabel('frequency (Hz)'); ylabel('phase coherence');

colorCode = colormap(jet(length(cellModels)));

plot(freq(:), meanSpect1(:,1))
figure;
hold on;
for i = 1:length(colorCode)
    plot(freq(:), meanSpect2(:,i), 'Color', colorCode(i,:))
end

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

plot(f, meanStimFFT, f, meanResidFFT)

infoSpect = log2(meanStimFFT./meanResidFFT);

plot(f, infoSpect)

% compute info spect density between stim and decoded residuals

Fs = 1000;
numWindows = 4;
nfft = 2^nextpow2(slen)/numWindows;
f = Fs/2*linspace(0,1,nfft/2+1);
% [signalSpect, freqs] = pwelch(outSignal, 4, [], [], 1000);

reshapedStimMat = reshape(outSignal,slen/numWindows,numSignals*numWindows);
StimFFT = fft(reshapedStimMat,nfft)/slen;
meanStimFFT = mean(2*abs(StimFFT(1:nfft/2+1,:)),2);

clear infoSpect
for i = 1:numPops
    for j = 1:numSignals
        reconStim(:,j) = reconStructOut(j+(i-1)*numSignals).optStim;
    end
    reconResids = reconStim-outSignal;
    reshapedResidMat = reshape(reconResids,slen/numWindows,numSignals*numWindows);
    residFFT = fft(reshapedResidMat,nfft)/slen;
    meanResidFFT = mean(2*abs(residFFT(1:nfft/2+1,:)),2);
    
    infoSpect(:,i) = log2(meanStimFFT./meanResidFFT);
end

% plot(f, meanStimFFT, f, meanResidFFT)

% infoSpect = log2(meanStimFFT./meanResidFFT);
allCellInds = [9]; %[10 18 33 40];
maxf = 27;

colorCode = colormap(jet(length(cellModels)));
figure;
hold on;
for i = 1:length(allCellInds)
    plot(f(1:maxf), infoSpect(1:maxf, allCellInds(i)), 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 2);
end
ylabel('info spectrum (bits/Hz)'); xlabel('frequency (Hz)');


figure; 
% allCellInds = unique(nonzeros(currCellInds));
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

%% get spectrum values for peak, width, etc
spectPeakVals = max(infoSpect',[], 2);
spectPeaks = f(argmax(infoSpect'));

% at what frequency does spect return to x pct of maximum?
returnPct = .3;

startFreq = zeros(1,numModels);
for i = 1:numModels
    startFreq(i) = f(find(infoSpect(1:end,i) > spectPeakVals(i)*returnPct, 1, 'first'));
    startFreqVal(i) = infoSpect(find(infoSpect(1:end,i) > spectPeakVals(i)*returnPct, 1, 'first'),i);
end

returnFreq = zeros(1,numModels);
for i = 1:numModels
%     validFreqs = f(infoSpect(:,i) >= startFreqVal(i));
%     validFreqs = f(2:end);
%     returnFreq(i) = f(find(infoSpect(3:end,i) < spectPeakVals(i)*returnPct, 1, 'first'));
    returnFreq(i) = f(find( (infoSpect(3:end,i) < spectPeakVals(i)*returnPct) & (f(3:end) > startFreq(i))', 1, 'first'));
    returnFreqVal(i) = infoSpect(find(infoSpect(3:end,i) < spectPeakVals(i)*returnPct & (f(3:end) > startFreq(i))', 1, 'first'),i);
end

spectWidth = abs(startFreq - returnFreq);
figure; plot(1:44, spectWidth,'.');

figure; hist(spectPeaks); 
ylabel('count'); xlabel('frequency (Hz)');

figure; hist(spectWidth, 5);
ylabel('count'); xlabel('frequency (Hz)');

allCellInds = find(spectWidth < 20);

allCellInds = [30]; %[10 18 33 40];
maxf = 27;

colorCode = colormap(jet(length(cellModels)));
figure;
hold on;
for i = 1:length(allCellInds)
    currCell = allCellInds(i);
    plot(f(1:maxf), infoSpect(1:maxf, allCellInds(i)), 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 2);
    plot([startFreq(currCell) returnFreq(currCell)], [startFreqVal(currCell) returnFreqVal(currCell)], '*');
end

%% check the infospectrum for pops of multiple cells
popIndMat = zeros(numModels, numModels);
cnt = 1;
for i = 1:numModels
    for j = i:numModels
        popIndMat(i,j) = cnt + numModels;
        popIndMat(j,i) = cnt + numModels;
        cnt = cnt + 1;
    end
end

% plot a histogram of diffEntVals for a given cell
currCell = 10;
compCell = 33;

currCellInds = [currCell 0; compCell 0; currCell currCell; compCell compCell; currCell compCell];

figure; hold on;
for i = 1:size(currCellInds,1)
    tempCurrCellInds = currCellInds(i,:);
    if nnz(tempCurrCellInds) < 2
        plotInds(i) = tempCurrCellInds(1);
        currColor = colorCode(tempCurrCellInds(1),:);
        currSymbol = '-';
%         symbol{i} = '-';
        plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), currSymbol, 'Color', currColor, 'LineWidth' , 2);
        
    elseif tempCurrCellInds(1) == tempCurrCellInds(2) %two cells and homo
        plotInds(i) = popIndMat(tempCurrCellInds(1), tempCurrCellInds(2));
        currColor = colorCode(tempCurrCellInds(1),:);
        currSymbol = '--';
        plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), currSymbol, 'Color', currColor, 'LineWidth' , 2);
%         symbol{i} = '--';
    else %two cells and hetero
        plotInds(i) = popIndMat(tempCurrCellInds(1), tempCurrCellInds(2));
        currColor = colorCode(tempCurrCellInds(1),:);
        currSymbol = '-';
        plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), currSymbol, 'Color', currColor, 'LineWidth' , 2.5);
        currColor = colorCode(tempCurrCellInds(2),:);
        currSymbol = '--';
        plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), currSymbol, 'Color', currColor, 'LineWidth' , 1.5);
%         symbol{i} = '--';
    end
%     colorCodePlot(i,:) = colorCode(tempCurrCellInds(1),:);
end 
ylabel('info spectrum (bits/Hz)');
% subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('frequency (Hz)');

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
% figure; hold on;
% for i = 1:size(currCellInds,1)
%     tempCurrCellInds = currCellInds(i,:);
%     plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), symbol{i}, 'Color', colorCodePlot(i,:), 'LineWidth' , 2)
% end

% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('info spectrum (bits/Hz)');
% subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('frequency (Hz)');

%% check the infospectrum for pops of multiple cells using data from
% pairwiseSims1

popIndMat = zeros(length(useCells), length(useCells));
cnt = 1;
for i = 1:length(useCells)
    for j = i:length(useCells)
        popIndMat(i,j) = cnt + length(useCells);
        popIndMat(j,i) = cnt + length(useCells);
        cnt = cnt + 1;
    end
end

% plot a histogram of diffEntVals for a given cell
currCell = 19;
compCell = 22;

currCellIndex = find(currCell == useCells);
compCellIndex = find(compCell == useCells);

currCellInds = [currCellIndex 0; compCellIndex 0; currCellIndex currCellIndex; compCellIndex compCellIndex; currCellIndex compCellIndex];
colorCellInds = [currCell 0; compCell 0; currCell currCell; compCell compCell; currCell compCell];


figure; hold on;
for i = 1:size(currCellInds,1)
    tempCurrCellInds = currCellInds(i,:);
    tempColorCellInds = colorCellInds(i,:);
    
    if nnz(tempCurrCellInds) < 2
        plotInds(i) = tempCurrCellInds(1);
        currColor = colorCode(tempColorCellInds(1),:);
        currSymbol = '-';
%         symbol{i} = '-';
        plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), currSymbol, 'Color', currColor, 'LineWidth' , 2);
        
    elseif tempCurrCellInds(1) == tempCurrCellInds(2) %two cells and homo
        plotInds(i) = popIndMat(tempCurrCellInds(1), tempCurrCellInds(2));
        currColor = colorCode(tempColorCellInds(1),:);
        currSymbol = '--';
        plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), currSymbol, 'Color', currColor, 'LineWidth' , 2);
%         symbol{i} = '--';
    else %two cells and hetero
        plotInds(i) = popIndMat(tempCurrCellInds(1), tempCurrCellInds(2));
        currColor = colorCode(tempColorCellInds(1),:);
        currSymbol = '-';
        plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), currSymbol, 'Color', currColor, 'LineWidth' , 2.5);
        currColor = colorCode(tempColorCellInds(2),:);
        currSymbol = '--';
        plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), currSymbol, 'Color', currColor, 'LineWidth' , 1.5);
%         symbol{i} = '--';
    end
%     colorCodePlot(i,:) = colorCode(tempCurrCellInds(1),:);
end 
ylabel('info spectrum (bits/Hz)');
% subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('frequency (Hz)');

figure; 
allCellInds = unique(nonzeros(colorCellInds));
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
% figure; hold on;
% for i = 1:size(currCellInds,1)
%     tempCurrCellInds = currCellInds(i,:);
%     plot(f(1:maxf), infoSpect(1:maxf, plotInds(i)), symbol{i}, 'Color', colorCodePlot(i,:), 'LineWidth' , 2)
% end

% % legend('real stim', 'same recon1', 'same recon2','diff recon');
% ylabel('info spectrum (bits/Hz)');
% % subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
% xlabel('frequency (Hz)');

%% are the max freq widths and heights for hetero > homo?

