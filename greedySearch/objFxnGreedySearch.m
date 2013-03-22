function [popAvgVal, popValsAll] = objFxnGreedySearch(reconStructOut, popAssignInds, stimParms, varargin)
numStimStats = stimParms.numStimStats;
numPops = max(popAssignInds.popInd);

if varargin{1} == 1 || 4
    switchInd = 1; % use mean recon error averaged over stimStats
    outSignal = varargin{2};
else
    switchInd = 2;
end

if switchInd==1
    meanReconAccAll = zeros(popAssignInds.totalInds, 1);
    parfor i = 1:popAssignInds.totalInds
        % first compute errors
        reconSignal = reconStructOut(i).optStim;
        %     reconSignal = reconStructOut(i).linStim;
        realSignal = outSignal(:,popAssignInds.stimInd(i));
        meanReconAccAll(i) = sqrt(mean((realSignal - reconSignal).^2));
        
    end
    sampleVals = meanReconAccAll;
else
    %stimEnt needs to get passed in
    %stimEnt = .5*sum(myLog2(eig(stimCovMat))) + .5*slen*myLog2(2*pi*exp(1));
    respEntSample = zeros(popAssignInds.totalInds, 1);
    parfor i = 1:popAssignInds.totalInds
        respEntSample(i) = -.5*sum(myLog2(reconStructOut(i).hessianEigs));
    end
%     totalEnt = stimEnt - (respEntSample  + .5*slen*myLog2(2*pi*exp(1)));
%     totalEntPerSec = totalEnt/(slen*.001);
%     sampleVals = totalEntPerSec;
    sampleVals = respEntSample;
end

popValsAll = [];
for j = 1:numPops
    tempValsAll = [];
    for k = 1:numStimStats
        %                 currInds = find(popAssignInds.divInd == i & popAssignInds.popInd == j);
        currInds = find(popAssignInds.popInd == j & popAssignInds.stimStatInd == k);
        popAvgVal(k,j) = mean(sampleVals(currInds));
        tempValsAll = vectCat(tempValsAll, sampleVals(currInds));
    end
    popValsAll = vectCat(popValsAll, tempValsAll);
end


%
% % for i = 1:length(divFact)
% for i = 1:length(corrVals)
%     for j = 1:numPops
% %                 currInds = find(popAssignInds.divInd == i & popAssignInds.popInd == j);
%         currInds = find(popAssignInds.corrInd == i & popAssignInds.popInd == j);
%         popMeanErrors(i,j) = mean(meanReconAccAll(currInds));
%     end
% end
%
%
% stimEnt = .5*sum(myLog2(eig(stimCovMat))) + .5*slen*myLog2(2*pi*exp(1));
% respEntSample = zeros(popAssignInds.totalInds, 1);
% parfor i = 1:popAssignInds.totalInds
%     respEntSample(i) = -.5*sum(myLog2(reconStructOut(i).hessianEigs));
% end
% totalEnt = stimEnt - (respEntSample  + .5*slen*myLog2(2*pi*exp(1)));
% totalEntPerSec = totalEnt/(slen*.001);
%
% % for i = 1:length(divFact)
% for i = 1:length(corrVals)
%     for j = 1:numPops
% %                 currInds = find(popAssignInds.divInd == i & popAssignInds.popInd == j);
%         currInds = find(popAssignInds.corrInd == i & popAssignInds.popInd == j);
%         popMeanEnt(i,j) = mean(totalEntPerSec(currInds));
%     end
% end

