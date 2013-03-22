
currPopSize = 1;
%find indexes of populations of currPopSize

cnt = 0;
for i = 1:length(pop)
    if length(pop(i).dat) == currPopSize
        relPopInds(cnt+1) = i;
        cnt = cnt + 1;
    end
end

%for each of the pops, compute their rmse's for each of the stims
% stimEnt = .5*sum(myLog2(eig(stimCovMat))) + .5*slen*myLog2(2*pi*exp(1));
for i = 1:length(relPopInds)
%     hessDetSample = zeros(1, numSignals);
    for j = 1:numSignals
        meanReconAccSample(j) = sqrt(mean((outSignal(:,j) - reconStructOut(j+(relPopInds(i)-1)*numSignals).optStim).^2));
    end
    reconRMSEAll(:,i) = meanReconAccSample;
%     reconRMSEMean(i) = mean(meanReconAccSample);
%     
    respEntSample = zeros(1, numSignals);
    for j = 1:numSignals
%         hessMatTemp = reconStructOut(j+(i-1)*numSignals).hessian;
        respEntSample(j) = -.5*sum(myLog2(reconStructOut(j+(relPopInds(i)-1)*numSignals).hessianEigs));
    end
    respEntAll(:,i) = respEntSample + .5*slen*myLog2(2*pi*exp(1));
% 
%     respEnt(i) = mean(respEntSample) + .5*slen*myLog2(2*pi*exp(1));
%     respEntError(i) = std(respEntSample)/sqrt(numSignals);
%     totalEnt(i) = stimEnt - respEnt(i);
%     totalEntError(i) = respEntError(i);
end

for i = 1:numStimStats
    stimEnt(i) = .5*sum(myLog2(eig(stimCovMat(:,:,i)))) + .5*slen*myLog2(2*pi*exp(1));
end

for i = 1:numStimStats
    reconRMSEMean(i,:) = mean(reconRMSEAll(1+(i-1)*numSignalsPerStat:i*numSignalsPerStat,:));
    respEntMean(i,:) = mean(respEntAll(1+(i-1)*numSignalsPerStat:i*numSignalsPerStat,:));
    totalEntMean(i,:) = stimEnt(i) - respEntMean(i,:);
end

totalEntPerSec = totalEntMean/(slen*.001);

