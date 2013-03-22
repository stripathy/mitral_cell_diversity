function [reconStructOut, decodingTime] = reconRunfile(cellModels,reconStructIn, stimCovMat, stimuli)

outSignal = stimuli;
%% get decoded responses
slen = size(outSignal,1);
stimInterval = [0 slen];
numSignals = size(outSignal,2);
% currStimCovMat = stimCovMat(:,:,1);
numSignalsPerStat = size(stimCovMat,3);
dtStim = 1;
tic
parfor i = 1:length(reconStructIn)
    
    if isfield(reconStructIn, 'stimInd');
        currStimInd = reconStructIn(i).stimInd;
    else
        currStimInd = mod(i-1,numSignals)+1;
    end
    initStim = outSignal(:,currStimInd);
    
    if isfield(reconStructIn, 'noisyInds'); 
        currPop = cellModels(reconStructIn(i).noisyInds);
    else
        currPop = cellModels(reconStructIn(i).popMakeup);
    end
    
    if numSignalsPerStat == 1 && ~isfield(reconStructIn, 'stimCovMatInd')
        currStimCovMat = stimCovMat(:,:,1);
    else
        covMatInd = reconStructIn(i).stimCovMatInd;
        if covMatInd ~= -1
            currStimCovMat = stimCovMat(:,:,covMatInd);
        else
            currStimCovMat = -1;
        end
    end
    
    testSpikes = reconStructIn(i).testSpikes;  
%     optFilts = reconStructIn(i).optFilts;
%     optStim = stimRecon(testSpikes, [1 slen], optFilts, staWin);
%     
    [optStim, exitflag, likeli, hessian] = bayesStimDecoder1(currPop, testSpikes, stimInterval, dtStim, currStimCovMat, initStim);
    reconStructOut(i).optStim = optStim;
    reconStructOut(i).likeli = likeli;
    reconStructOut(i).hessianEigs = eig(hessian);
    reconStructOut(i).stimErrorBars = sqrt(diag(inv(hessian)));
%     reconStructOut(i).exitflag = exitflag;
    if currStimCovMat == -1
        reconStructOut(i).hessian = sparse(triu(hessian,0));
    end
end
decodingTime = toc/60;

