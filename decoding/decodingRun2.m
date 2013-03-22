%% example of decoding done on simulated spikes from 3 pops:
% 2 pops of homo cells and 1 hetero pop

cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP');
load('LNP_recons_drugs5.mat');
load('alphaStimCovMat');

pop(1).cell(1:5) = ggAll(1);
pop(2).cell(1:5) = ggAll(13);
pop(3).cell(1:5) = ggAll([1 13 15 12 8]);

inputStim = getNoisyStim(250, 1, 3);
stimInterval = [1 250];
numRepeats = 50;
optStim = zeros(250, numRepeats*length(pop));
clear reconAccVec testSpikesAll
parfor k = 1:numRepeats*length(pop)
    currTrialInd = mod(k-1,numRepeats) + 1;
    if k < numRepeats + 1
        j = 1;
    elseif k < 2*numRepeats + 1
        j = 2;
    else
        j = 3;
    end
   
        testSpikes = [];
            for i = 1:length(pop(j).cell)
                [tempSpikes] = simGLM(pop(j).cell(i), inputStim);
                tempSpikes = tempSpikes(tempSpikes>stimInterval(1) & tempSpikes<stimInterval(2));
                tempSpikes = tempSpikes -stimInterval(1) + 1;
                testSpikes = vectCat(testSpikes, tempSpikes);
            end
            testSpikesAll(k).dat = testSpikes;
            if (isempty(testSpikes))
                continue;
            end
            
            cellModels = pop(j).cell;
            spikeTimes = testSpikes;
            dtStim = 1; %ms
        
        % stimCovMat = eye(stimInterval(2)-stimInterval(1) + 1);
        initStim = rand(length(inputStim), 1);
        [optStim(:,k), exitflag(k), fval(k)] = bayesStimDecoder1(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat, initStim);
%         reconAccTemp = reconAcc(inputStim, optStim(:,k));
        
%     reconAccVec(k) = reconAccTemp;
end

for i = 1:numRepeats*length(pop)
    [reconAccVec(i) reconInfo(i)] = reconAcc(inputStim, optStim(:,i));
end

reconAccAll = reshape(reconAccVec', numRepeats, length(reconAccVec)/numRepeats);
reconInfo = reshape(reconInfo', numRepeats, length(reconInfo)/numRepeats);

groupErrors = mean(reconAccAll);
groupStdErr = std(reconAccAll)./sqrt(numRepeats);
figure;barweb(groupErrors, groupStdErr);

groupErrors = mean(reconInfo);
groupStdErr = std(reconInfo)./sqrt(numRepeats);
figure;barweb(groupErrors, groupStdErr); 

[h12 p12] = ranksum(reconAccAll(:,1), reconAccAll(:,2));
[h23 p23] = ranksum(reconAccAll(:,2), reconAccAll(:,3));
[h13 p13] = ranksum(reconAccAll(:,1), reconAccAll(:,3));

%plot an example trace with spikes
colorInds = ['bgr'];
rasterColorInds(1:5) = colorInds(1);
rasterColorInds(6:10) = colorInds(2);
rasterColorInds(11:15) = colorInds(3);

figure; subplot(4,1,1:2); plot(1:250, inputStim, 'k', 'LineWidth' , 2)
axis([1 250 -3 3]);
currTrialInd = 5;
plotInds = numRepeats*[0:2] + currTrialInd;
testSpikesPlot = [];
hold on;
for i =1:3
    subplot(4,1,1:2); plot(1:250, optStim(:,plotInds(i)), colorInds(i), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, testSpikesAll(plotInds(i)).dat);
    subplot(4,1,3); hold on; plot(1:250, inputStim-optStim(:,plotInds(i)), colorInds(i), 'LineWidth' , 2)
    axis([1 250 -2.5 2.5]);
end
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('stimulus value (zscore)');
subplot(4,1,4);createRaster2(testSpikesPlot',1, 250, rasterColorInds); ylabel('Neuron ID');
xlabel('time (ms)');

tptsK = -length(sta)+1:1:0;
clear staMat kMat ihMat
figure
colorInds = ['bgr'];
for i = 1:length(pop)
    for j = 1:length(pop(i).cell)
        staMat(:,j) = pop(i).cell(j).sta;
        kMat(:,j) = pop(i).cell(j).k;
        ihMat(:,j) = pop(i).cell(j).ihbas*pop(i).cell(j).ih;
    end
    
%     subplot(3,2,2*i-2);plot(1:length(sta), staMat, colorInds(i)); axis([0 50 -.5 1.5]);
    subplot(3,2,2*i-1);plot(tptsK, kMat, colorInds(i) , 'LineWidth' , 2); axis([-50 0 -.5 1.5]);
    subplot(3,2,2*i);plot(ggAll(1).iht,ihMat, colorInds(i), 'LineWidth' , 2); axis([0 63 -2 3.2]);
end