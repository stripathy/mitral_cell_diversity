stimIndex = repmat(1:10, 5, 1);
stimIndex = stimIndex(:);
stimIndex = repmat(stimIndex, 2, 1);


cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\greedySearch\finSims\sim3');

% load data from greedySearch iterations across stimuli
% these get stored in bestPop

fileLists = dir('greedySearch*');
[fileLists] = sortGreedySearchNames(fileLists);

nameBase = 'greedySearch1';
typeBase = '.mat';

numStimStats = length(fileLists);
popCnt = 1;
popListBest = zeros(5, numStimStats*10);
bestPop = zeros(10, 100);
for s = 1:numStimStats
    load(fileLists(s).name);
    temp = strrep(fileLists(s).name, nameBase, []);
    num = str2num(strrep(temp, typeBase, []));
%     load(['greedySearch1',num2str(s)]);
    tempBestPop = getFinalPop(bestPopListAll);
    [bestPop(:,num)] = tempBestPop;
    
%     for i = 1:length(tempBestPop)
%         popListBest(1:i, popCnt) = sort(tempBestPop(1:i));
%         popCnt = popCnt + 1;
%     end
end

cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\greedySearch\finSims\sim5');

% load data from greedySearch iterations across stimuli
% these get stored in bestPop

fileLists = dir('greedySearch*');
[fileLists] = sortGreedySearchNames(fileLists);

nameBase = 'greedySearch1';
typeBase = '.mat';

numStimStats = length(fileLists);
popCnt = 1;
popListBest = zeros(5, numStimStats*10);
% bestPop = [];
for s = 1:numStimStats
    load(fileLists(s).name);
    temp = strrep(fileLists(s).name, nameBase, []);
    num = str2num(strrep(temp, typeBase, []))+50;
%     load(['greedySearch1',num2str(s)]);
    tempBestPop = getFinalPop(bestPopListAll);
    [bestPop(:,num)] = tempBestPop;
    
%     for i = 1:length(tempBestPop)
%         popListBest(1:i, popCnt) = sort(tempBestPop(1:i));
%         popCnt = popCnt + 1;
%     end
end

numStimStats = length(fileLists);

numBestPops = size(bestPop,2);
clear *Diff;
for i = 1:numBestPops
    tempPopList = nonzeros(bestPop(:,i));
    
    [stimDiff(:,i) histDiff(:,i) biasDiff(:,i)] = glmDiff(cellModels(tempPopList));
%     [stimDiffUni(:,i) histDiffUni(:,i) biasDiffUni(:,i)] = glmDiff(cellModels(unique(tempPopList)));
end

stimTypes = 1:10;
clear *DiffAll
for i = 1:length(stimTypes)
    currStimType = stimTypes(i);
    inds = find(stimIndex == i);
    stimDiffAll(:, i) = stimDiff(inds);
    histDiffAll(:, i) = histDiff(inds);
    biasDiffAll(:, i) = biasDiff(inds);
end

bestAvgDivStim = nanmean(stimDiffAll(:));
bestAvgDivHist = nanmean(histDiffAll(:));
bestAvgDivBias = nanmean(biasDiffAll(:));

stdErrDivStim = myStdErr(stimDiffAll(~isnan(stimDiffAll)));
stdErrDivHist = myStdErr(histDiffAll(~isnan(histDiffAll)));
stdErrDivBias = myStdErr(biasDiffAll(~isnan(biasDiffAll)));

stdErrDivStim = myStdErr(stimDiffAll(1,:));
stdErrDivHist = myStdErr(histDiffAll(1,:));
stdErrDivBias = myStdErr(biasDiffAll(1,:));


% stimDiffAll = reshape(stimDiff,numSearchIters,numStimStats)';
% histDiffAll = reshape(histDiff,numSearchIters,numStimStats)';
% biasDiffAll = reshape(biasDiff,numSearchIters,numStimStats)';

removeStims = [7 9];
includeStims = setdiff(1:10, removeStims);
numStimStats = length(fileLists);
figure;
for k = 1:3
    if k == 1
        useMat = stimDiffAll;
    elseif k == 2
        useMat = histDiffAll;
    else
        useMat = biasDiffAll;
    end
    p = zeros(10,10);
    for i = 1:10
        vec1 = useMat(:, i);
        vec1 = vec1(~isnan(vec1));
        for j = i+1:10
            vec2 = useMat(:, j);
            vec2 = vec2(~isnan(vec2));
            p(i,j) = ranksum(vec1, vec2);
        end
    end
p = p + p';

subplot(1,3,k);
plotMat = p(includeStims, includeStims);
imagesc(log10(plotMat)); colorbar; caxis([-2 0]);
hold on;
[r, c] = find(plotMat < .05 & plotMat ~= 0);
plot(c, r, 'ow')
% figure;
% bar(mean(useMat));
% 
% hold on;
% errorbar(1:10, mean(useMat), myStdErr(useMat), 'r.');
end

        
plotStims = [4 5 2];
figure; hold on;
for k = 1:3
    if k == 1
        useMat = stimDiffAll;
        meanDiv = meanDivStim(10);
    elseif k == 2
        useMat = histDiffAll;
        meanDiv = meanDivHist(10);
    else
        useMat = biasDiffAll;
        meanDiv = meanDivBias(10);
    end

    xPts = [ones(10,1).*randn(10,1)*.1 ones(10,1).*randn(10,1)*.1 + 1 ones(10,1).*randn(10,1)*.1 + 2];
    xPts = xPts + (k-1)*4;
%     plot(xPts, useMat(:,plotStims)/meanDiv, '.');
    plot(xPts(:,1), useMat(:,plotStims(1))/meanDiv, 'm.');
    plot(xPts(:,2), useMat(:,plotStims(2))/meanDiv, 'c.');
    plot(xPts(:,3), useMat(:,plotStims(3))/meanDiv, 'k.');
    
end
axis([-1 11 -.2 1.5])

%% 
hetErrorBars = [(meanDivStim(10) - meanCIStim(10,1)), (meanDivHist(10) - meanCIHist(10,1)),...
    (meanDivBias(10) - meanCIBias(10,1))]./[meanDivStim(10), meanDivHist(10), meanDivBias(10)];
bestErrorBars = [stdErrDivStim, stdErrDivHist,stdErrDivBias]./[meanDivStim(10), meanDivHist(10), meanDivBias(10)];
vals = [bestAvgDivStim, bestAvgDivHist, bestAvgDivBias]./[meanDivStim(10), meanDivHist(10), meanDivBias(10)];
figure;
plot([1 2 3], [1 1 1], 'g.', [1 2 3], [0 0 0], 'r.', [1 2 3], vals, 'b.');
hold on;
errorbar([1 2 3], [1 1 1], hetErrorBars, '.g')
errorbar([1 2 3], vals, bestErrorBars, '.b')
axis([0 4 -.2 1.5])

