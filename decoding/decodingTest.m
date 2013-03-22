cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP');
load('LNP_recons_3_05_10.mat');

stimInterval = [2301 2550];

tempStim = getNoisyStim(100000, 1, 3);
crossCorr = xcov(tempStim, 50, 'coeff'); %
crossCorr = crossCorr(51:101);
crossCorr(15:end) = 0;
% crossCorr(crossCorr<0) = 0;
crossCorr(end+1:stimInterval(2)-stimInterval(1) + 1) = 0;
stimCovMat = toeplitz(crossCorr);

%start with data for one cell
cellInds = 1;
currTrialInd = 3;

testSpikes = [];
for i = 1:length(cellInds)
    tempSpikes = spikeTimeMat(:,ggAllCells(cellInds(i)).trialInds(currTrialInd));
    tempSpikes = tempSpikes(tempSpikes>stimInterval(1) & tempSpikes<stimInterval(2));
    testSpikes = vectCat(testSpikes, tempSpikes);
end

cellModels = ggAllCells(cellInds);
spikeTimes = testSpikes;
dtStim = 1; %ms

% stimCovMat = eye(stimInterval(2)-stimInterval(1) + 1);

[optStim, exitflag, logli, hessian] = bayesStimDecoder(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat);

invHess = inv(hessian);
errorBarsOptStim = sqrt(diag(invHess));

figure; errorbar(1:250, optStim, errorBarsOptStim);

[eigVects eigSpect] = eig(invHess);
plot(sqrt(diag(eigSpect)), '.')

[minVal, minInd] = min(errorBarsOptStim);
plot(eigVects(:,2))

tPts = stimInterval(1):dtStim:stimInterval(2);
spikeVec = histc(spikeTimes, tPts);
figure; plot(1:250, testStim, 'b', 1:250, optStim, 'g', 1:250, spikeVec-1, 'r')

plot(1:50, crossCorr(1:50)); xlabel('time'); ylabel('correlation value');
% figure; plot(1:250, testStim, 1:250, optStim(:,exitflag==2))
%% run on LNP neuron

slen = 250; % Stimulus length (frames) 
swid = 1;  % Stimulus width  (pixels).  Must match # pixels in stim filter
dtStim = 1;
testStim = randn(slen,swid);  % Gaussian white noise stimulus
[testSpikes] = simGLM(ggsim, testStim);  % Simulate GLM response
stimInterval = [1 250];
stimCovMat = eye(stimInterval(2)-stimInterval(1) + 1);

[optStim] = bayesStimDecoder(gg, testSpikes, stimInterval, dtStim, stimCovMat);

% figure; plot(1:250, testStim, 1:250, optStim)
figure; plot(1:250, testStim, 1:250, optStim, 'r')
%% replicate linear decoding figure using Bayesian
cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP');
load('LNP_recons_10bases.mat');

stimInterval = [2301 2550];

groupInds = [1 4 3];
currTrialInd = 18;

% cellInds = [1 7];
for j = 1:length(groupInds)
    for k = 1:length(allGroupInds(j,:))
        cellInds(j,k) = find(allGroupInds(j,k) == cellId);
        %         cellInds(k) = find(allGroupInds(1,k) == cellId);
    end
end

testSpikesAll = [];
for j = 1:length(groupInds)
%     for k = 1:length(allGroupInds(j,:))
%         cellInds(k) = find(allGroupInds(j,k) == cellId);
% %         cellInds(k) = find(allGroupInds(1,k) == cellId);
%     end
    
    testSpikes = [];
    for i = 1:length(cellInds)
        tempSpikes = spikeTimeMat(:,cellStartInds(cellId(cellInds(j,i)))+currTrialInd - 1);
        tempSpikes = tempSpikes(tempSpikes>stimInterval(1) & tempSpikes<stimInterval(2));
        tempSpikes = tempSpikes -stimInterval(1) + 1;
        testSpikes = vectCat(testSpikes, tempSpikes);
    end
    
    cellModels = ggAll(cellInds);
    spikeTimes = testSpikes;
    dtStim = 1; %ms
    testSpikesAll = vectCat(testSpikesAll, testSpikes);

    % stimCovMat = eye(stimInterval(2)-stimInterval(1) + 1);
    initStim = rand(length(testStim), 1);
    [optStim(:,j), exitflag(j), fval(j)] = bayesStimDecoder1(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat, initStim);
end

figure; plot(1:250, testStim, 1:250, optStim(:,exitflag~=0), 'r')

colorInds = ['rgk'];
rasterColorInds(1:5) = colorInds(1);
rasterColorInds(6:10) = colorInds(2);
rasterColorInds(11:15) = colorInds(3);
figure; subplot(3,1,1:2); plot(2301:2300+250, testStim)
hold on;
for i =1:3
    plot(2301:2300+250, optStim(:,i), colorInds(i))
end
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus value (zscore)');
subplot(3,1,3);createRaster2(testSpikesAll',1, 250, rasterColorInds);
xlabel('time (ms)');


figure
groupInds = [1 4 3];
currCellIds = [];
colorInds = ['rgk'];
for i = 1:length(groupInds)
    for j = 1:5
        currCellIds(j) = find(allGroupInds(i,j)==cellId);
        staMat(:,j) = ggAll(currCellIds(j)).sta;
        kMat(:,j) = ggAll(currCellIds(j)).k;
        ihMat(:,j) = ggAll(currCellIds(j)).ihbas*ggAll(currCellIds(j)).ih;
    end
    
    subplot(3,3,3*i-2);plot(1:length(sta), staMat, colorInds(i)); axis([0 50 -.5 1]);
    subplot(3,3,3*i-1);plot(1:length(sta), kMat, colorInds(i)); axis([0 50 -.5 1]);
    subplot(3,3,3*i);plot(ggAll(1).iht,ihMat, colorInds(i)); axis([0 70 -10 4]);
end

%% testing multiple runs of the same neuron(s)
cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP');
load('LNP_recons_10bases.mat');
load('alphaStimCovMat');

stimInterval = [2301 2450];
dtStim = 1;
tPts = stimInterval(1):dtStim:stimInterval(2);
slen = length(tPts);

tempStim = getNoisyStim(2000000, dtStim, 3);
% crossCorr = xcov(tempStim, staWin/dtStim, 'coeff'); %make sure this is monotonically decreasing
crossCorr = xcov(tempStim, slen, 'coeff');
crossCorr = crossCorr(slen+1:slen*2);
% plot(1:101, crossCorr, 1:.5:101, crossCorr2)

crossCorr = crossCorr(staWin/dtStim + 1: 2*staWin/dtStim + 1);
crossCorr = crossCorr(staWin/dtStim + 1: 2*staWin/dtStim + 1);
firstNegInd = find(diff(crossCorr)>0, 1);
crossCorr(firstNegInd:end) = 0;
crossCorr(crossCorr<0) = 0;

crossCorr = exp(-.5*(0:dtStim:50));

crossCorr(end+1:slen) = 0;
stimCovMat = toeplitz(crossCorr);


stimCovMat = eye(slen);

[mat, corrmat] = corrmtx(tempStim, slen-1);

% tauAlpha = 3; dt = .001; tptsAlpha = dt:dt:15;
% filt = tptsAlpha.*exp(-tptsAlpha/tauAlpha)/tauAlpha;
% y = fft(filt);
% ccorr = ifft(abs(y).^2);

% groupInds = [1];
currTrialInd = 18;
numRepeats = 1;

cellInds(1) = 5;
inputStim = randn(slen, 1);
% inputStim = getNoisyStim(slen*dtStim,1, 3);
% [testSpikes] = simGLM(ggAll(cellInds), inputStim) + stimInterval(1);
testSpikes = [];
for i = 1:length(cellInds)
    [tempSpikes] = simGLM(ggAll(cellInds(i)), inputStim);
    testSpikes = vectCat(testSpikes, tempSpikes);
end

tempStim = randn(slen,1);

% tempStim = getNoisyStim(slen*dtStim,dtStim, 3);

%get the right cellInds for the group
% for j = 1:length(groupInds)
%     for k = 1:length(allGroupInds(j,:))
%         cellInds(j,k) = find(allGroupInds(j,k) == cellId);
%         %         cellInds(k) = find(allGroupInds(1,k) == cellId);
%     end
% end

testSpikesAll = [];
optStim = [];
exitflag = []; fval = [];

for j = 1:numRepeats
%     for k = 1:length(allGroupInds(j,:))
%         cellInds(k) = find(allGroupInds(j,k) == cellId);
% %         cellInds(k) = find(allGroupInds(1,k) == cellId);
%     end
    
%     testSpikes = [];
%     for i = 1:length(cellInds)
%         tempSpikes = spikeTimeMat(:,cellStartInds(cellId(cellInds(i)))+currTrialInd - 1);
% %         tempSpikes = spikeTimeMat(:,cellStartInds(cellId(cellInds(j,i)))+currTrialInd - 1);
%         tempSpikes = tempSpikes(tempSpikes>stimInterval(1) & tempSpikes<stimInterval(2));
%         testSpikes = vectCat(testSpikes, tempSpikes);
%     end

% generate testSpikes from model itself
    
    
    cellModels = ggAll(cellInds);
    spikeTimes = testSpikes;
%     dtStim = .1; %ms
    testSpikesAll = vectCat(testSpikesAll, testSpikes);

    % stimCovMat = eye(stimInterval(2)-stimInterval(1) + 1);
    
    [optStim(:,j), exitflag(j), fval(j)] = bayesStimDecoder1(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat, tempStim);
end

tPts = stimInterval(1):dtStim:stimInterval(2);
spikeVec = histc(spikeTimes, tPts);
figure; plot(stimInterval(1):1:stimInterval(2), testStim, 'b', tPts, optStim, 'g')

figure
plot(tPts(1):1:tPts(end), inputStim, tPts, optStim, tPts, histc(spikeTimes, tPts))

%check that my fxn and pillow neglogli give similar values
[optStim(:,j), exitflag(j), fval(j)] = bayesStimDecoder1(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat, tempStim);
ggTest = cellModels(1);
ggTest.tsp = nonzeros(testSpikes);
neglogli_GLM(ggTest,tempStim)

%% working on getting bayesStimDecoderLogli non-prior part to give unique
% solns

dtStim = 1;
msInput = 250;
stimInterval = [1 msInput];
stimCovMat = eye(msInput,msInput);

currTrialInd = 18;
numRepeats = 50;

cellInds = [5 5 5 5 5 5 5 6 9 8 7 10 15 16 1 2 3];
% inputStim = randn(msInput, 1);
inputStim = getNoisyStim(msInput,dtStim, 3);
testSpikes = [];
for i = 1:length(cellInds)
    [tempSpikes] = simGLM(ggAll(cellInds(i)), inputStim);
    testSpikes = vectCat(testSpikes, tempSpikes);
end

% tempStim = randn(msInput, 1);
tempStim = fastinterp2(inputStim,1/dtStim);

tempStim = zeros(msInput*1/dtStim, 1);
cnt = 1;
for i = 1:length(inputStim)
    tempStim(cnt:cnt+1/dtStim - 1) = inputStim(i);
    cnt = cnt+1/dtStim;
end

% tempStim = getNoisyStim(slen*dtStim,dtStim, 3);
% tempStim = inputStim;
tempStim = randn(msInput,1);

testSpikesAll = [];
optStim = [];
exitflag = []; fval = [];

parfor j = 1:numRepeats
%     tempStim = randn(msInput,1);
    cellModels = ggAll(cellInds);
    spikeTimes = testSpikes;
    
    [optStim(:,j), exitflag(j), fval(j), hessian] = bayesStimDecoder1(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat, tempStim);
end

plot(1:msInput, inputStim, 1:msInput, optStim)

% test neglogli
ggTest = cellModels(1);
ggTest.tsp = nonzeros(testSpikes);
[neglogli, rr, tt, Iinj,Istm,Ih] = neglogli_GLM(ggTest,inputStim);

figure; plot(tt, Istm, tt, Ih)

%% 
invHess = inv(hessian);
[eVects eVals] = eig(invHess);
plot(diag(eVals),'.')
plot(eVects(:,2))

ebars = sqrt(diag(eVals));
errorbar(optStim(:,1), ebars)