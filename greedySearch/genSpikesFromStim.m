function [simSpikes, simSpikesInds] = genSpikesFromStim(outSignal, stimParms, varargin)

numStimStats = stimParms.numStimStats;
% varargin = cellModels, numStimStats

if isempty(varargin)
    % load cellModels
    load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\drugAllFitsFinalggAll.mat';
    cellModels = ggAll;
    numStimStats = 1;
    keepCells = 1:length(ggAll);
else
    cellModels = varargin{1};
    keepCells = 1:length(cellModels);
%     numStimStats = varargin{2};
end

slen = size(outSignal,1);
dtStim = 1;
numSignals = size(outSignal,2);
numSignalsPerStat = numSignals/numStimStats;


numSimRepeats = 10;

testInterval = [0 slen*dtStim];
maxSpikes = 100;
if slen > 512
    maxSpikes = 300;
end

% set up a list of indexes for simulating spike trains
% list things to iterate through:
% numSignals numModels divFact

clear simSpikesInds
simSpikesInds.stimInd = zeros(1, numSignals*length(keepCells));
simSpikesInds.modelInd = zeros(1, numSignals*length(keepCells));
simSpikesInds.stimStatInd = zeros(1, numSignals*length(keepCells));
cnt = 1;
for i = 1:numSignals
    for j = 1:length(keepCells) %make sure keepCells matches what's in cellModels
        simSpikesInds.stimInd(cnt) = i;
        simSpikesInds.modelInd(cnt) = keepCells(j); % this corresponds to the abstract model index
        simSpikesInds.stimStatInd(cnt) = (floor( (i-1)/numSignalsPerStat)+1);
        cnt = cnt + 1;
    end
end
simSpikesInds.numSimRepeats = numSimRepeats;

clear simSpikes
parfor i = 1:length(simSpikesInds.stimInd)
    currStimInd = simSpikesInds.stimInd(i);
    currModelInd = simSpikesInds.modelInd(i);
    
    currStim = outSignal(:,currStimInd);
    estSpikes = zeros(maxSpikes, numSimRepeats); %init spike matrix
    
    
    for j = 1:numSimRepeats
        [tspEst] = simGLM(cellModels(currModelInd), currStim); % cellModelsAct here is the "real" neuron model
        %         [tspEst] = simGLM(cellModelsAct(currModelInd), currStim);
        estSpikes(1:length(tspEst),j) = tspEst;
    end
    currSpikes = getSpikesInt(estSpikes, testInterval); % get spiketrains in correct time interval
    if currSpikes==0
        currSpikes = zeros(1, numSimRepeats);
    end
    [simSpikes(i).dat] = currSpikes;
end
