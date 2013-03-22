function [popAvgVal, popValsAll, popStdVal] = computeAvgErrors(reconStructOut, popAssignInds, stimParms, outSignal)
numStimStats = stimParms.numStimStats;
numPops = max(popAssignInds.popInd);

switchInd = 1;
% if varargin{1} == 1
%     switchInd = 1; % use mean recon error averaged over stimStats
%     outSignal = varargin{2};
% else
%     switchInd = 2; % use info (bits/sec)
% %     stimEnt = varargin{2};
% %     slen = size(stimParms.outSignal,1);
% end

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
popStdVal = [];
for j = 1:numPops
    tempValsAll = [];
    for k = 1:numStimStats
        %                 currInds = find(popAssignInds.divInd == i & popAssignInds.popInd == j);
        currInds = find(popAssignInds.popInd == j & popAssignInds.stimStatInd == k);
        popAvgVal(k,j) = mean(sampleVals(currInds));
        popStdVal(k,j) = std(sampleVals(currInds));
        tempValsAll = vectCat(tempValsAll', sampleVals(currInds)')';
    end
    popValsAll = vectCat(popValsAll, tempValsAll);
end

popAvgVal = popAvgVal';
popStdVal = popStdVal';
