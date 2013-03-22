% combines template1 with stuff wtih diff stats

numSignalsOld = numSignals;
outSignalOld = outSignal;
stimCovMatOld = stimCovMat;
simSpikesOld = simSpikes;
numStimStatsOld = numStimStats;

numStimStats = numStimStatsOld + numStimStats;
numSignals = numSignalsOld + numSignals;

stimCovMatTemp = stimCovMat;
stimCovMat = stimCovMatOld;
stimCovMat(:,:,4) = stimCovMatTemp;

simSpikesTemp = simSpikes;
simSpikes = simSpikesOld;
simSpikes(end+1:end+length(simSpikesTemp)) = simSpikesTemp;

outSignalTemp = outSignal;
outSignal = vectCat(outSignalOld, outSignalTemp);

save 'template2'
