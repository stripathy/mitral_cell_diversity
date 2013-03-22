% popInfoSim to simulate populations of varying sizes and using stimuli
% with different statistics

%% generate a bunch of different stims
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\diffStimStats\stimCovMats\finStimCovMat4';
crossCorrAll = crossCorrAll([1],:);
% generate stims
slen = 1024;
numSignalsPerStat = 1;
randSeedFixed = 0;

[outSignal, stimCovMat, stimParms] = generateStims(crossCorrAll, slen, numSignalsPerStat, randSeedFixed);

% generate spike trains
% load neuron models
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\drugAllFitsFinalggAll.mat';
keepCells = 1:length(ggAll);
numModels = length(keepCells);
for i = 1:numModels
    cellModels(i) = ggAll(keepCells(i));
end

[simSpikes, simSpikesInds] = genSpikesFromStim(outSignal, stimParms, cellModels);
firingRateNeuronComp

% make populations to search
neuronsToSim = [1:44]; %don't use 1 or 7
maxPopSize = 2;

clear pop;
cnt = 1;
for i = 1:length(neuronsToSim)
    for j = maxPopSize
        pop(cnt).dat = repmat(neuronsToSim(i), 1, j);
        cnt = cnt + 1;
    end
end

numPopsPerSize = 1;
for i = 1:maxPopSize
    for j = 1:numPopsPerSize
        pop(cnt).dat = sort(randsample(keepCells(neuronsToSim), i ,'true')');
        cnt = cnt + 1;
    end
end

numPops = length(pop);

% code to assign spikes to pops and make reconStructIn

[reconStructIn, popAssignInds] = assignPopSpikes(pop, simSpikes, simSpikesInds, stimParms);

% 
    

