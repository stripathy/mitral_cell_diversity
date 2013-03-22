% merge data from cluster split


numSplits = 10;
cnt = 1;
for i = 1:numSplits
    loadFile = ['split' num2str(i) 'Out.mat'];
    if exist(loadFile,'file') == 0
        continue;
    end
    load(loadFile);
    inds = (cnt-1)*length(reconStructOut) + 1:cnt*length(reconStructOut);
%     reconStructInTemp(inds) = currReconStructIn;
    reconStructOutTemp(inds) = reconStructOut;
    cnt = cnt + 1;
end

clear reconStructOut %reconStructIn;
% unrandomize indexes
% [blah, randInds] = sort(rand(length(pop),1));
for i = 1:length(randInds)
    structInInds = (i-1)*numSignals+1:i*numSignals;
    structOutInds = (randInds(i)-1)*numSignals+1:randInds(i)*numSignals;
    reconStructOut(structOutInds) = reconStructOutTemp(structInInds);
%     reconStructIn(structOutInds) = reconStructInTemp(structInInds);
end

% reconStructOut = reconStructOutTemp;
clear reconStructTemp currReconStructIn reconStructOutTemp currInds reconStructInTemp

%% merge data across multiple files

folderDir = 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\diffStimStats\';

folders{1} = 'sim2';
folders{2} = 'sim3';
folders{3} = 'sim4';
folders{4} = 'sim7';

cd ([folderDir folders{1}]);
load 'reconsMerged';

neuronsToSimTemp = [];

numPopsTemp = numPops;
popTemp = pop;
reconStructInTemp = reconStructIn;
reconStructOutTemp = reconStructOut;
neuronsToSimTemp = union(neuronsToSim, neuronsToSimTemp);

for i = 2:length(folders)
    cd ([folderDir folders{i}]);
    load 'reconsMerged';
    
    numPopsTemp = numPopsTemp + numPops;
    popTemp(length(popTemp)+1:length(popTemp)+length(pop)) = pop;
    
    reconStructInTemp(length(reconStructInTemp)+1:length(reconStructInTemp)+length(reconStructIn)) = reconStructIn;
    reconStructOutTemp(length(reconStructOutTemp)+1:length(reconStructOutTemp)+length(reconStructOut)) = reconStructOut;
    neuronsToSimTemp = union(neuronsToSim, neuronsToSimTemp);
end
    
numPops = numPopsTemp;
pop = popTemp;
reconStructIn = reconStructInTemp;
reconStructOut = reconStructOutTemp;
neuronsToSim = neuronsToSimTemp;

clear numPopsTemp popTemp reconStructInTemp reconStructOutTemp neuronsToSimTemp;

cd (folderDir);
save '5NeuronPopsWRecons';

clear reconStructIn reconStructOut
save '5NeuronPopsNoRecons';