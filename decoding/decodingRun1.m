% decoding run 1

%% test of accuracy of Bayesian Decoding vs Linear decoding for populations

cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP');
load('LNP_recons_10bases.mat');
load('alphaStimCovMat');

stimInterval = [2301 2550];
dtStim = 1;
tPts = stimInterval(1):dtStim:stimInterval(2);
slen = length(tPts);

groupInds = [1 4 3];
currTrialInd = 15;

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
    
    cellModels = ggAll(cellInds(j,:));
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
    plot(2301:2300+250, optStim(:,i), colorInds(i), 'LineWidth' , 2)
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

% for each of the groups and each trial, run the decoding on actual spikes
stimInterval = [2301 2550];
parfor k = 1:40*length(groupInds)
    currTrialInd = mod(k-1,40) + 1;
    if k < 41
        j = 1;
    elseif k < 81
        j = 2;
    else
        j = 3;
    end
   
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
        
        cellModels = ggAll(cellInds(j,:));
        spikeTimes = testSpikes;
        dtStim = 1; %ms
        
        testSpikesAll(k).dat = testSpikes;
        
        % stimCovMat = eye(stimInterval(2)-stimInterval(1) + 1);
        initStim = rand(length(testStim), 1);
        [optStim(:,k), exitflag(k), fval(k)] = bayesStimDecoder1(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat, initStim);
        reconAccTemp = reconAcc(testStim, optStim(:,k));
        
    reconAccVec(k) = reconAccTemp;
end 

%plot an example trace with spikes
colorInds = ['rgk'];
rasterColorInds(1:5) = colorInds(1);
rasterColorInds(6:10) = colorInds(2);
rasterColorInds(11:15) = colorInds(3);
figure; subplot(3,1,1:2); plot(2301:2300+250, testStim, 'LineWidth' , 2)
currTrialInd = 15;
plotInds = numRepeats*[0:2] + currTrialInd;
testSpikesPlot = [];
hold on;
for i =1:3
    plot(2301:2300+250, optStim(:,plotInds(i)), colorInds(i), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, testSpikesAll(plotInds(i)).dat);
end
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus value (zscore)');
subplot(3,1,3);createRaster2(testSpikesPlot',1, 250, rasterColorInds);
xlabel('time (ms)');

groupErrors = mean([reconAccVec(1:40); reconAccVec(41:80); reconAccVec(81:120)]');
groupStdErr = std([reconAccVec(1:40); reconAccVec(41:80); reconAccVec(81:120)]')./sqrt(40);
figure;barweb(groupErrors, groupStdErr);

%% run decoding on simualated spikes generated from the model


% for each of the groups and each trial, run the decoding on actual spikes
% inputStim = testStim;
inputStim = getNoisyStim(250, 1, 3);
stimInterval = [1 250];
numRepeats = 100;
optStim = zeros(250, numRepeats*length(groupInds));
clear reconAccVec testSpikesAll
parfor k = 1:numRepeats*length(groupInds)
    currTrialInd = mod(k-1,numRepeats) + 1;
    if k < numRepeats + 1
        j = 1;
    elseif k < 2*numRepeats + 1
        j = 2;
    else
        j = 3;
    end
   
        testSpikes = [];
        inputCellIndVec = [];
            for i = 1:length(cellInds(j,:))
                inputCellInd = i;
                if j == 1
                    inputCellInd = 1;
                elseif j == 2
                    inputCellInd = 3;
                else
                    inputCellInd = i;
                end
                inputCellIndVec(i) = inputCellInd;
                [tempSpikes] = simGLM(ggAll(cellInds(j,inputCellInd)), inputStim);
                tempSpikes = tempSpikes(tempSpikes>stimInterval(1) & tempSpikes<stimInterval(2));
                tempSpikes = tempSpikes -stimInterval(1) + 1;
                testSpikes = vectCat(testSpikes, tempSpikes);
            end
            testSpikesAll(k).dat = testSpikes;
            if (isempty(testSpikes))
                continue;
            end
            
            cellModels = ggAll(cellInds(j, inputCellIndVec));
            spikeTimes = testSpikes;
            dtStim = 1; %ms
        
        % stimCovMat = eye(stimInterval(2)-stimInterval(1) + 1);
        initStim = rand(length(testStim), 1);
        [optStim(:,k), exitflag(k), fval(k)] = bayesStimDecoder1(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat, initStim);
        reconAccTemp = reconAcc(inputStim, optStim(:,k));
        
    reconAccVec(k) = reconAccTemp;
end

reconAccAll = reshape(reconAccVec', numRepeats, length(reconAccVec)/numRepeats);

groupErrors = mean(reconAccAll);
groupStdErr = std(reconAccAll)./sqrt(numRepeats);
figure;barweb(groupErrors, groupStdErr);

[h12 p12] = ranksum(reconAccAll(:,1), reconAccAll(:,2));
[h23 p23] = ranksum(reconAccAll(:,2), reconAccAll(:,3));
[h13 p13] = ranksum(reconAccAll(:,1), reconAccAll(:,3));

%plot an example trace with spikes
colorInds = ['rgk'];
rasterColorInds(1:5) = colorInds(1);
rasterColorInds(6:10) = colorInds(2);
rasterColorInds(11:15) = colorInds(3);
figure; subplot(3,1,1:2); plot(1:250, inputStim, 'LineWidth' , 2)
currTrialInd = 8;
plotInds = numRepeats*[0:2] + currTrialInd;
testSpikesPlot = [];
hold on;
for i =1:3
    plot(1:250, optStim(:,plotInds(i)), colorInds(i), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, testSpikesAll(plotInds(i)).dat);
end
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus value (zscore)');
subplot(3,1,3);createRaster2(testSpikesPlot',1, 250, rasterColorInds);
xlabel('time (ms)');
