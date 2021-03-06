% popInfoSim to simulate populations of varying sizes

%% generate a bunch of different stims
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\decoding\alphaStimCovMat';

slen = 256;
slenBurn = 50;
corrVec = stimCovMat(1,1:14);
corrVec(end+1:slen+slenBurn) = 0;
stimCovMat = toeplitz(corrVec);

% stimCovMat = eye(slen+slenBurn);
numSignals= 50;
dtStim = 1;

% cnt = 1;
% outSignal = zeros(slen, numSignals);
% for j = 1:numSignals
%     outSignal(:,j) = getNoisyStim(slen, 1, 3);
% end

cholMat = chol(stimCovMat);
outSignal = zeros(slen+slenBurn, numSignals);
for j = 1:numSignals
    outSignal(:,j) = cholMat*randn(slen+slenBurn,1);
end

outSignal = outSignal(slenBurn+1:end,:);
stimCovMat= stimCovMat(slenBurn+1:end, slenBurn+1:end);

% dtStim = 1;
if dtStim ~= 1
    tempOutSignal = zeros(slen*dtStim, numSignals);
    for i = 1:numSignals
        tempOutSignal(:,i) = reshape(repmat(outSignal(:,i),1,dtStim)', slen*dtStim, 1);
    end
    outSignal = tempOutSignal;
end



% clear corrVec
% numCorrPts = slen;
% tau = 1;
% corrVec(1:numCorrPts) = exp((-tau*(0:numCorrPts-1))).^5;

tapers=3; avg=1; fs = 1000;
pars = struct ('Fs', fs, ...
    'tapers', [tapers 2*tapers-1], ...
    'pad', 0, ...
    'trialave', 1);
for k = 1:numSignals
    [S(:,k), f] = mtspectrumc(outSignal(:,k), pars);
end
plot(f,mean(S'))


%% for every cell in my pop, stimulate each with every stim with x repeats

% load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\LNP_recons_drugs5.mat';
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\drugAllFitsFinalggAll.mat';

% keepCells = [1 2 7 11 8 14 3 13 12 6 15 5];
% keepCells = fliplr(keepCells);
keepCells = 1:44;

% [sortVals sortInds] = sort(meanFR);
% keepCells = sortInds;


numModels = length(keepCells);

for i = 1:numModels
    cellModels(i) = ggAll(keepCells(i));
end

numSimRepeats = 100;
testInterval = [0 slen*dtStim];
maxSpikes = 100;
if slen > 512
    maxSpikes = 300;
end


tic
clear simRates
parfor j = 1:numSignals*numModels
    modelInd = floor((j-1)/numSignals) +1;
    estSpikes = zeros(maxSpikes, numSimRepeats);
    for i = 1:numSimRepeats
        [tspEst] = simGLM(cellModels(modelInd), outSignal(:,mod(j-1,numSignals)+1));
        estSpikes(1:length(tspEst),i) = tspEst;
    end
    [simSpikes(j).dat] = getSpikesInt(estSpikes, testInterval);
    
end
time1 = toc;

%% construct pops from cell-spike responses and decode pop spike trains
% neuronsToSim = [7 18 32 35 41];
% neuronsToSim = [11:21];
% neuronsToSim = [27:35];
neuronsToSim = [38:43];
maxPopSize = 40;

clear pop;
cnt = 1;
for i = 1:length(neuronsToSim)
    for j = 1:maxPopSize
        pop(cnt).dat = repmat(neuronsToSim(i), 1, j);
        cnt = cnt + 1;
    end
end

numPopsPerSize = 50;
for i = 1:maxPopSize
    for j = 1:numPopsPerSize
        pop(cnt).dat = sort(randsample(keepCells(neuronsToSim), i ,'true')');
        cnt = cnt + 1;
    end
end
% clear pop reconStructIn
% pop(1).dat = sort(repmat([1:44], 1, 150));
% pop(1).dat = 44;

numPops = length(pop);


for i = 1:length(pop)
    
    currCells = pop(i).dat;
    popModels = cellModels(currCells);
    for j = 1:numSignals
        tempTestSpikes = [];
        [uniqueCells, firstInds] = unique(currCells,'first');
        [uniqueCells, lastInds] = unique(currCells,'last');
        for k = 1:length(uniqueCells)
            cellStimInd = j + (uniqueCells(k)-1)*numSignals;
            randTrials = randsample(numSimRepeats, lastInds(k) - firstInds(k) + 1);
%             randTrials = 1;%randsample(numSimRepeats, length(currCells(k)));
%             if length(currCells) == 2
%                 if currCells(1) == currCells(2) && k == 2
%                     randTrials = 2;
%                 end
%             end
            tempCurrSpikes = simSpikes(cellStimInd).dat(:,randTrials);
            tempTestSpikes = vectCat(tempTestSpikes, tempCurrSpikes);
        end
        
        currPopInd = j + (i-1)*numSignals;
        reconStructIn(currPopInd).testSpikes = tempTestSpikes;
        reconStructIn(currPopInd).id = i;
        reconStructIn(currPopInd).popMakeup = currCells;
    end
end


%% get decoded responses
slen = size(outSignal,1);
stimInterval = [0 slen];
% dtStim = 2;
tic
parfor i = 1:length(reconStructIn)
    initStim = outSignal(:,mod(i-1,numSignals)+1);
    currPop = cellModels(reconStructIn(i).popMakeup);
    testSpikes = reconStructIn(i).testSpikes;  
%     optFilts = reconStructIn(i).optFilts;
%     optStim = stimRecon(testSpikes, [1 slen], optFilts, staWin);
%     
    [optStim, exitflag, likeli, hessian] = bayesStimDecoder1(currPop, testSpikes, stimInterval, dtStim, stimCovMat, initStim);
    reconStructOut(i).optStim = optStim;
    reconStructOut(i).likeli = likeli;
    reconStructOut(i).hessianEigs = eig(hessian);
    reconStructOut(i).stimErrorBars = sqrt(diag(inv(hessian)));
%     reconStructOut(i).exitflag = exitflag;
    reconStructOut(i).hessian = hessian;
end
decodingTime = toc/60;

figure;
stairs(1:length(optStim), downsample(outSignal,dtStim), 'k'); 
hold on;
stairs(1:length(optStim), optStim, 'r');