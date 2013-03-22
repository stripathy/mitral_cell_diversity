%% generate a bunch of different stims
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\decoding\alphaStimCovMat';

slen = 1024;
slenBurn = 50;
corrVec = stimCovMat(1,1:14);
corrVec(end+1:slen+slenBurn) = 0;
stimCovMat = toeplitz(corrVec);

numSignals= 400;

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
figure; plot(f,mean(S')); ylabel('Power spectra'); xlabel('Frequency (Hz)');
axis([0 100 0 .012]); box off;



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
testInterval = [0 slen];
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

% useCells = [3 8 10 18 19 22 31 33 42 43];
useCells = [keepCells];

% clear pop;
% for i = 1:numModels
%     pop(i).dat = [i];
% end

% cnt = length(pop)+1;
% for i = 1:numModels
%     for j = i:numModels
%         pop(cnt).dat = [i, j];
%         cnt = cnt + 1;
%     end
% end

clear pop;
for i = 1:length(useCells)
    pop(i).dat = useCells(i);
end

cnt = length(pop)+1;
for i = 1:length(useCells)
    for j = i:length(useCells)
        pop(cnt).dat = [useCells(i), useCells(j)];
        cnt = cnt + 1;
    end
end

numPops = length(pop);


for i = 1:length(pop)
    
    currCells = pop(i).dat;
    popModels = cellModels(currCells);
    for j = 1:numSignals
        tempTestSpikes = [];
        for k = 1:length(currCells)
            cellStimInd = j + (currCells(k)-1)*numSignals;
            %             randTrials = randsample(numSimRepeats, length(currCells(k)));
            randTrials = 1;%randsample(numSimRepeats, length(currCells(k)));
            if length(currCells) == 2
                if currCells(1) == currCells(2) && k == 2
                    randTrials = 2;
                end
            end
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
dtStim = 1;
tic
parfor i = 1:length(reconStructIn)
    initStim = outSignal(:,mod(i-1,numSignals)+1);
    currPop = cellModels(reconStructIn(i).popMakeup);
    testSpikes = reconStructIn(i).testSpikes;  
%     optFilts = reconStructIn(i).optFilts;
%     optStim = stimRecon(testSpikes, [1 slen], optFilts, staWin);
%     
    [optStim, exitflag, likeli, hessian] = bayesStimDecoder1(currPop, testSpikes, stimInterval, dtStim, stimCovMat(:,:,1), initStim);
    reconStructOut(i).optStim = optStim;
    reconStructOut(i).likeli = likeli;
    reconStructOut(i).hessianEigs = eig(hessian);
    reconStructOut(i).stimErrorBars = sqrt(diag(inv(hessian)));
%     reconStructOut(i).exitflag = exitflag;
%     reconStructOut(i).hessian = hessian;
end
decodingTime = toc/60;