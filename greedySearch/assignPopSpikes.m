function [reconStructIn, popAssignInds] = assignPopSpikes(popStruct, simSpikes, simSpikesInds, stimParms)

numSignals = stimParms.numSignals;
numStimStats = stimParms.numStimStats;
% numSimRepeats = simSpikesInds.numSimRepeats;
numSignalsPerStat = numSignals/numStimStats;

numPops = length(popStruct);

clear popAssignInds;
cnt = 1;
popAssignInds.totalInds = numPops*numSignals;
popAssignInds.stimInd = zeros(1, popAssignInds.totalInds);
popAssignInds.popInd = zeros(1, popAssignInds.totalInds);
popAssignInds.stimStatInd = zeros(1, popAssignInds.totalInds);


for j = 1:numPops
    for k = 1:numSignals
        popAssignInds.stimInd(cnt) = k;
        popAssignInds.stimStatInd(cnt) = (floor( (k-1)/numSignalsPerStat)+1);
        popAssignInds.popInd(cnt) = j;
        %             popAssignInds.neuronIds(cnt).dat = (i-1)*length(keepCells)+ pop(j).dat;
        cnt = cnt + 1;
    end
end
% popAssignInds.totalInds = cnt-1;

% assign vals to reconStructIn: spike trains and indices to which cell models
clear reconStructIn
parfor i = 1:popAssignInds.totalInds
    currNeuronIds = popStruct(popAssignInds.popInd(i)).dat;
    %     currNeuronIds = popAssignInds.neuronIds(i).dat;
    currStimInd = popAssignInds.stimInd(i);
    
    [uniqueCells, firstInds] = unique(currNeuronIds,'first');
    [uniqueCells, lastInds] = unique(currNeuronIds,'last');
    tempTestSpikes = [];
    for k = 1:length(uniqueCells)
        currSpikesInd = find(uniqueCells(k) == simSpikesInds.modelInd & currStimInd == simSpikesInds.stimInd);
        
        if simSpikesInds.numSimRepeats == -1 % then each neuron has an unequal number of repeats
            numSimRepeats =  size(simSpikes(currSpikesInd).dat,2);
        else
            numSimRepeats = simSpikesInds.numSimRepeats;
        end
        randTrials = randsample(numSimRepeats, lastInds(k) - firstInds(k) + 1);
        tempCurrSpikes = simSpikes(currSpikesInd).dat(:,randTrials);
        tempTestSpikes = vectCat(tempTestSpikes, tempCurrSpikes);
    end
    
    reconStructIn(i).testSpikes = tempTestSpikes;
    reconStructIn(i).id = i;
    reconStructIn(i).popMakeup = currNeuronIds;
    reconStructIn(i).stimCovMatInd = (floor( (currStimInd-1)/numSignalsPerStat)+1);
    reconStructIn(i).stimInd = currStimInd;
end

% 
% numPops = length(popStruct);
% for i = 1:numPops
%     
%     currCells = popStruct(i).dat;
%     popModels = cellModels(currCells);
%     for j = 1:numSignals
%         tempTestSpikes = [];
%         [uniqueCells, firstInds] = unique(currCells,'first');
%         [uniqueCells, lastInds] = unique(currCells,'last');
%         for k = 1:length(uniqueCells)
%             cellStimInd = j + (uniqueCells(k)-1)*numSignals;
%             randTrials = randsample(numSimRepeats, lastInds(k) - firstInds(k) + 1); %sample spikeVecs w/o replacement
%             %             randTrials = 1;%randsample(numSimRepeats, length(currCells(k)));
%             %             if length(currCells) == 2
%             %                 if currCells(1) == currCells(2) && k == 2
%             %                     randTrials = 2;
%             %                 end
%             %             end
%             tempCurrSpikes = simSpikes(cellStimInd).dat(:,randTrials);
%             tempTestSpikes = vectCat(tempTestSpikes, tempCurrSpikes);
%         end
%         
%         currPopInd = j + (i-1)*numSignals;
%         reconStructIn(currPopInd).testSpikes = tempTestSpikes;
%         reconStructIn(currPopInd).id = i;
%         reconStructIn(currPopInd).popMakeup = currCells;
%         reconStructIn(currPopInd).stimCovMatInd = (floor( (j-1)/numSignalsPerStat)+1);
%     end
% end