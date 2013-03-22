% popInfoSim to simulate populations of varying sizes and using stimuli
% with different statistics

%% generate a bunch of different stims
% load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\alphaStimCovMatHighLow';
% load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\natStimCovMat2';
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\oscStimCovMat';
% crossCorrAll = crossCorrAll(2,:);

slen = 512;
slenBurn = 100;
numSignalsPerStat = 50;
numStimStats = size(crossCorrAll,1);

numSignals= numSignalsPerStat*numStimStats;
dtStim = 1;


stimCovMat = zeros(slen+slenBurn, slen+slenBurn, numStimStats);
outSignal = zeros(slen+slenBurn, numSignals);
for k = 1:numStimStats
    corrVec = nonzeros(crossCorrAll(k,:))';
    
    if length(corrVec) < slen+slenBurn
        corrVec(end+1:slen+slenBurn) = 0;
    else
        corrVec = corrVec(1:slen+slenBurn);
    end
    stimCovMat(:,:,k) = toeplitz(corrVec);
    
    
    cholMat = chol(stimCovMat(:,:,k));
    
    % not supporting dt's that aren't 1 ms yet
    for j = 1:numSignalsPerStat
        outSignal(:,j+(k-1)*numSignalsPerStat) = cholMat*randn(slen+slenBurn,1);
    end
    
end
outSignal = outSignal(slenBurn+1:end,:);
stimCovMat= stimCovMat(slenBurn+1:end, slenBurn+1:end,:);

stimParms.outSignal = outSignal;
stimParms.numSignals = numSignals;
stimParms.numStimStats = numStimStats;

for i = 1:numStimStats
    stimEnt(i) = .5*sum(myLog2(eig(stimCovMat(:,:,i)))) + .5*slen*myLog2(2*pi*exp(1));
end
% 
tapers=3; avg=1; fs = 1000;
pars = struct ('Fs', fs, ...
    'tapers', [tapers 2*tapers-1], ...
    'pad', 0, ...
    'trialave', 1);
for k = 1:numStimStats
    [S(:,k), f] = mtspectrumc(outSignal(:,1+(k-1)*numSignalsPerStat:k*numSignalsPerStat), pars);
end
figure; plot(f,S); ylabel('power'); xlabel('Frequency (Hz)');
axis([0 100 0 .04]);
loglog(f,S); ylabel('power'); xlabel('Frequency (Hz)');
plot(f,S)

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

numSimRepeats = 10;
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
%         if isempty(tspEst);
%             tspEst = 0;
%         end
        estSpikes(1:length(tspEst),i) = tspEst;
    end
    currSpikes = getSpikesInt(estSpikes, testInterval);
     if isequal(currSpikes, 0)
        currSpikes = zeros(2, numSimRepeats);
    end   
    
    [simSpikes(j).dat] = currSpikes;
    
end
time1 = toc;

%% construct pops from cell-spike responses and decode pop spike trains
% neuronsToSim = [7 18 32 35 41];
% neuronsToSim = [11:21];
% neuronsToSim = [27:35];
% neuronsToSim = [38:43];
neuronsToSim = [2:6 8:44]; %don't use 1 or 7

% pop(1).dat = [30 33 41 43 44];
% pop(2).dat = [27 31 33 36 44];
% pop(3).dat = [27 31 33 36 44];
% pop(4).dat = [27 31 33 41 44];

maxPopSize = 5;

clear pop;
cnt = 1;
for i = 1:length(neuronsToSim)
    for j = maxPopSize
        pop(cnt).dat = repmat(neuronsToSim(i), 1, j);
        cnt = cnt + 1;
    end
end


numPopsPerSize = 4000;
% tempPopMat = zeros(maxPopSize, numPopsPerSize);
for i = maxPopSize
    for j = 1:numPopsPerSize
        pop(cnt).dat = sort(randsample(keepCells(neuronsToSim), i ,'true')');
%         tempPopMat(:,j) = pop(cnt).dat;
        cnt = cnt + 1;
    end
end

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
%             if isequal(simSpikes(cellStimInd).dat,0)
%                 tempCurrSpikes = zeros(length(randTrials),1);
%                 
%             else
            tempCurrSpikes = simSpikes(cellStimInd).dat(:,randTrials);
%             end
            tempTestSpikes = vectCat(tempTestSpikes, tempCurrSpikes);
        end
        
        currPopInd = j + (i-1)*numSignals;
        reconStructIn(currPopInd).testSpikes = tempTestSpikes;
        reconStructIn(currPopInd).id = i;
        reconStructIn(currPopInd).popMakeup = currCells;
        reconStructIn(currPopInd).stimCovMatInd = (floor( (j-1)/numSignalsPerStat)+1);
    end
end


%% get decoded responses
slen = size(outSignal,1);
stimInterval = [0 slen];
% currStimCovMat = stimCovMat(:,:,1);
numSignalsPerStat = size(stimCovMat,3);
% dtStim = 2;
tic
parfor i = 1:length(reconStructIn)
    initStim = outSignal(:,mod(i-1,numSignals)+1);
    currPop = cellModels(reconStructIn(i).popMakeup);
    
    if numSignalsPerStat == 1
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
%     reconStructOut(i).likeli = likeli;
%     reconStructOut(i).hessianEigs = eig(hessian);
%     reconStructOut(i).stimErrorBars = sqrt(diag(inv(hessian)));
% %     reconStructOut(i).exitflag = exitflag;
%     reconStructOut(i).hessian = sparse(hessian);
end
decodingTime = toc/60;

% figure;
% stairs(1:length(optStim), downsample(outSignal,dtStim), 'k'); 
% hold on;
% stairs(1:length(optStim), optStim, 'r');