%% generate a bunch of different stims
load 'alphaStimCovMat';

slen = 256;
slenBurn = 50;
corrVec = stimCovMat(1,1:14);
corrVec(end+1:slen+slenBurn) = 0;
stimCovMat = toeplitz(corrVec);

numStims= 1;

% cnt = 1;
% stim = zeros(slen, numStims);
% for j = 1:numStims
%     stim(:,j) = getNoisyStim(slen, 1, 3);
% end

cholMat = chol(stimCovMat);
stim = zeros(slen+slenBurn, numStims);
for j = 1:numStims
    stim(:,j) = cholMat*randn(slen+slenBurn,1);
end

stim = stim(slenBurn+1:end,:);
stimCovMat= stimCovMat(slenBurn+1:end, slenBurn+1:end);

% clear corrVec
% numCorrPts = slen;
% tau = 1;
% corrVec(1:numCorrPts) = exp((-tau*(0:numCorrPts-1))).^5;

% plot powerspectrum of inputs if you have chronux toolbox
% tapers=3; avg=1; fs = 1000;
% pars = struct ('Fs', fs, ...
%     'tapers', [tapers 2*tapers-1], ...
%     'pad', 0, ...
%     'trialave', 1);
% for k = 1:numStims
%     [S(:,k), f] = mtspectrumc(stim(:,k), pars);
% end
% plot(f,mean(S'))


%% for every cell in my pop, stimulate each with every stim with x repeats

% load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\LNP_recons_drugs5.mat';
load 'drugAllFitsFinalggAll';

% keepCells = [1 2 7 11 8 14 3 13 12 6 15 5];
% keepCells = fliplr(keepCells);
keepCells = 1:44;

% [sortVals sortInds] = sort(meanFR);
% keepCells = sortInds;


numModels = length(keepCells);

for i = 1:numModels
    cellModels(i) = ggAll(keepCells(i));
end

numSimRepeats = 50;
testInterval = [0 slen];
maxSpikes = 100;


tic
clear simRates
parfor j = 1:numStims*numModels
    modelInd = floor((j-1)/numStims) +1;
    estSpikes = zeros(maxSpikes, numSimRepeats);
    for i = 1:numSimRepeats
        [tspEst] = simGLM(cellModels(modelInd), stim(:,mod(j-1,numStims)+1));
        estSpikes(1:length(tspEst),i) = tspEst;
    end
    [simSpikes(j).dat] = getSpikesInt(estSpikes, testInterval);
    
end
time1 = toc;

%% construct pops from cell-spike responses and decode pop spike trains

%example pop

popMakeup = [20 20]; %sort(randi(44,5,1));
pop.dat = popMakeup;

% for i = 1:numModels
%     pop(length(pop)+1).dat = [i, i];
% end
% for i = 1:numModels
%     pop(length(pop)+1).dat = [i, i, i, i, i];
% end

numPops = length(pop);


for i = 1:length(pop)
    
    currCells = pop(i).dat;
    popModels = cellModels(currCells);
    for j = 1:numStims
        tempTestSpikes = [];
        for k = 1:length(currCells)
            cellStimInd = j + (currCells(k)-1)*numStims;
            %             randTrials = randsample(numSimRepeats, length(currCells(k)));
            randTrials = randsample(numSimRepeats, length(currCells(k)));
            if length(currCells) == 2
                if currCells(1) == currCells(2) && k == 2
                    randTrials = 2;
                end
            end
            tempCurrSpikes = simSpikes(cellStimInd).dat(:,randTrials);
            tempTestSpikes = vectCat(tempTestSpikes, tempCurrSpikes);
        end
        
        currPopInd = j + (i-1)*numStims;
        reconStructIn(currPopInd).testSpikes = tempTestSpikes;
        reconStructIn(currPopInd).id = i;
        reconStructIn(currPopInd).popMakeup = currCells;
    end
end


%% get decoded responses
stimInterval = [0 slen];
dtStim = 1;
tic
parfor i = 1:length(reconStructIn)
    initStim = stim(:,mod(i-1,numStims)+1);
    currPop = cellModels(reconStructIn(i).popMakeup);
    testSpikes = reconStructIn(i).testSpikes;  
%     optFilts = reconStructIn(i).optFilts;
%     optStim = stimRecon(testSpikes, [1 slen], optFilts, staWin);
%     
    [optStim, exitflag, likeli, hessian] = bayesStimDecoder1(currPop, testSpikes, stimInterval, dtStim, stimCovMat, initStim);
    reconStructOut(i).optStim = optStim;
    reconStructOut(i).likeli = likeli;
    reconStructOut(i).hessianEigs = eig(hessian);
%     reconStructOut(i).exitflag = exitflag;
%     reconStructOut(i).hessian = hessian;
end
decodingTime = toc/60;


%% plot some recons and spike trains
colorCodeTotal = colormap(jet(length(keepCells)));

for i = 1:length(popMakeup)
    rasterColors(i,:) = colorCodeTotal(popMakeup(i),:);
end

% currCellInd = find(keepCells == myCellInd);
% currCellInd = 3;
currStimInd = randi(numStims,1);

figure; subplot(3,1,1:2); plot(1:slen, stim(:,currStimInd), 'k', 'LineWidth' , 2)
% axis([1 slen -3 3]);

plotInds = currStimInd;
testSpikesPlot = [];
hold on;
for i =1
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim, 'r', 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
end
axis tight;
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus (a.u.)'); set(gca,'XTick', []);
subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); ylabel('neuron');
xlabel('time (ms)');   

