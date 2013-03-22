cnt = 1;
for i = 1:stimParms.numStimStats
    for j = 1:numSignalsPerStat
        stimStatInds(cnt) = i;
        cnt = cnt + 1;
    end
end


currStimStat = 3;

stimInds = stimStatInds==currStimStat;
nfft = slen*1;
Fs = 1000;
f = Fs/2*linspace(0,1,nfft/2+1);


stim = outSignal(:,stimInds);
StimFFT = fft(stim,nfft)/slen;
meanStimFFT = mean(2*abs(StimFFT(1:nfft/2+1,:)),2);

meanReconFFT = zeros(length(f),numPops);
infoSpect = zeros(length(f),numPops);
for i = 1:numPops
    currPop = i;
    reconInds = find(popAssignInds.popInd == currPop & popAssignInds.stimStatInd== currStimStat);
    for j = 1:length(reconInds)
        reconStim(:,j) = reconStructOut(reconInds(j)).optStim;
    end
    reconResids = reconStim-stim;
    residFFT = fft(reconResids,nfft)/slen;
    meanResidFFT = mean(2*abs(residFFT(1:nfft/2+1,:)),2);
    
    reconFFT = fft(reconStim,nfft)/slen;
    meanReconFFT(:,i) = mean(2*abs(reconFFT(1:nfft/2+1,:)),2);
    
    infoSpect(:,i) = log2(meanStimFFT./meanResidFFT);
end

freqInds = 1:53;
close all
plot(f(freqInds), meanStimFFT(freqInds, :), 'k')
hold on; plot(f(freqInds), meanReconFFT(freqInds, :))
axis tight

freqPlot = 1:106;

figure;
plot(f(freqPlot), meanStimFFT(freqPlot), 'k', f(freqPlot), meanReconFFT(freqPlot,[44 155]));


% compute 95% conf intervals for the PSDs using residFFT
numResamps = 1000;
psdAll = zeros(length(1:nfft/2+1), numResamps);
for k = 1:numResamps
    resampInds = randi(numSignalsPerStat, 1, numSignalsPerStat);
    reconFFTResamp = reconFFT(:,resampInds );
    
    psdAll(:,k) = mean(2*abs(reconFFTResamp(1:nfft/2+1,:)),2);
end

for i = 1:size(psdAll,1)
    confInt(:,i) = getConfInterval(psdAll(i,:), 95);
end