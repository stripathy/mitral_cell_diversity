%% set up template files for greedy search

% go through diff stim cov mat files, initialize stims, load cellModels and
% which cells to use
stimDir = 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\diffStimStats\stimCovMats';
cd(stimDir);

stimFiles = {'alphaStimCovMatMultiple', 'ouStimCovMat', 'oscStimCovMat', 'natStimCovMat1', 'whiteStimCovMat'...
    'ouStimCovMat3', 'natStimCovMat5'};

stimFiles = {'finStimCovMat4'};
numStimFiles = length(stimFiles);


saveDir = 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\greedySearch\science_test\sim6';

keepStims = [1 2 3 4 5 6 7 9 10];
slen = 256;
numSignals = 20;
randSeedFixed = 0;
numStatRepeats = 1;


%load cellModels
load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\model fitting\drugAllFitsFinalggAll.mat';
keepCells = 1:length(ggAll);
numModels = length(ggAll);
for i = 1:numModels
    cellModels(i) = ggAll(keepCells(i));
end
neuronsToSim = [1:length(ggAll)]; %don't use 1 or 7

outSignal = [];
templateCnt = 1;
for i = 1:numStimFiles
    cd(stimDir)
    load(stimFiles{i});
    numStimStats = size(crossCorrAll,1);
    
    for j = 1:numStimStats
        for k = 1:numStatRepeats
            assert(randSeedFixed==0)
            %         load
            [outSignal, stimCovMat, stimParms] = generateStims(crossCorrAll(j,:), slen, numSignals, randSeedFixed);
            
%                     [outSignalTemp, stimCovMat, stimParms] = generateStims(crossCorrAll(j,:), slen, numSignals, randSeedFixed);
%                     outSignal = vectCat(outSignal, outSignalTemp);
            
            cd(saveDir);
            save(['template', num2str(templateCnt)]);
            templateCnt = templateCnt + 1;
        end
    end
end

% thing to do for multiple stims
outSignal = [];
templateCnt = 1;
for i = 1:numStimFiles
    cd(stimDir)
    load(stimFiles{i});
%     crossCorrAll = crossCorrAll([1:6 8 10],:);
    crossCorrAll = crossCorrAll(keepStims,:);
    
    numStimStats = size(crossCorrAll,1);
    
    assert(randSeedFixed==0)
    
    [outSignal, stimCovMat, stimParms] = generateStims(crossCorrAll, slen, numSignals, randSeedFixed);
    cd(saveDir);
%     save(['template', num2str(templateCnt)]);
    save(['template2']);
    templateCnt = templateCnt + 1;

end
        