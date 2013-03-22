%% generate a bunch of different stims
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\decoding\alphaStimCovMat';

slen = 3000; %stim length (tpts)
slenBurn = 50;
corrVec = stimCovMat(1,1:14);
corrVec(end+1:slen+slenBurn) = 0;
stimCovMat = toeplitz(corrVec);

numStims= 1;

cholMat = chol(stimCovMat);
stim = zeros(slen+slenBurn, numStims);
for j = 1:numStims
    stim(:,j) = cholMat*randn(slen+slenBurn,1);
end

stim = stim(slenBurn+1:end,:);
stimCovMat= stimCovMat(slenBurn+1:end, slenBurn+1:end);

%% plot the power spectrum of the stimulus if you have the chronux toolbox
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

% load the glm parameters for each neuron in the data set
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\drugAllFitsFinalggAll.mat';

keepCells = 1:44;

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
    currSpikes = getSpikesInt(estSpikes, testInterval);
    if isequal(currSpikes, 0)
        currSpikes = zeros(2, numSimRepeats);
    end   
    [simSpikes(j).dat] = currSpikes;
end
time1 = toc;

%% construct pops from cell-spike responses and decode pop spike trains

%example pop

popMakeup = [5 9 13 22 30];
pop.dat = popMakeup;

numPops = length(pop);

currCells = pop.dat;
popModels = cellModels(currCells);
for j = 1:numStims
    tempTestSpikes = [];
    for k = 1:length(currCells)
        cellStimInd = j + (currCells(k)-1)*numStims;
        randTrials = randsample(numSimRepeats, length(currCells(k)));
        tempCurrSpikes = simSpikes(cellStimInd).dat(:,randTrials);
        tempTestSpikes = vectCat(tempTestSpikes, tempCurrSpikes);
    end
    
    currPopInd = j + (i-1)*numStims;
    reconStructIn(currPopInd).testSpikes = tempTestSpikes;
    reconStructIn(currPopInd).id = i;
    reconStructIn(currPopInd).popMakeup = currCells;
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

%% compute mutual info - don' worry about this
% compute entropy using laplace approx on posterior MAP distn
stimEnt = .5*sum(myLog2(eig(stimCovMat))) + .5*slen*myLog2(2*pi*exp(1));
for i = 1:numPops
%     hessDetSample = zeros(1, numStims);
    respEntSample = zeros(1, numStims);
    for j = 1:numStims
%         hessMatTemp = reconStructOut(j+(i-1)*numStims).hessian;
        respEntSample(j) = -.5*sum(myLog2(reconStructOut(j+(i-1)*numStims).hessianEigs)) + .5*slen*myLog2(2*pi*exp(1));
    end
%     respEnt(i) = mean(-.5*myLog2(hessDetSample)) + .5*slen*myLog2(2*pi*exp(1));
    respEnt(i) = mean(respEntSample);
end

totalEnt = (stimEnt - respEnt);
totalEntPerSec = totalEnt/(slen*.001);

%% plot some recons and spike trains
colorCodeTotal = colormap(jet(length(keepCells)));

for i = 1:length(popMakeup)
    rasterColors(i,:) = colorCodeTotal(popMakeup(i),:);
end

currCellInd = find(keepCells == myCellInd);
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

%% plot a bar graph of stim vs population mutual info. - don't worry about this
stimEntPerSec = stimEnt/(slen*.001);

ciData = (stimEnt - respEntSample)/(slen*.001);
[totalEntCI] = myBootstrapCI(ciData, .95);
effErrorBars = abs(totalEntPerSec - totalEntCI);


figure;
bar(1, stimEntPerSec, 'FaceColor', 'k'); hold on
bar(2, totalEntPerSec, 'FaceColor', 'r');
set(gca,'XTick', []);
ylabel('information (bits/s)');
errorbar(2, totalEntPerSec, effErrorBars(1), effErrorBars(2), '-k', 'LineWidth', 1);

% figure
% barweb([stimEntPerSec totalEntPerSec], [);
