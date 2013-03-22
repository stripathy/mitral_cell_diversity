%% load the folder containing data
cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\greedySearch\finSims\sim2');

% load data from greedySearch iterations across stimuli
% these get stored in bestPop

fileLists = dir('greedySearch*');
[fileLists] = sortGreedySearchNames(fileLists);

numStimStats = length(fileLists);
popCnt = 1;
popListBest = zeros(5, numStimStats*10);
for s = 1:numStimStats
    load(fileLists(s).name);
%     load(['greedySearch1',num2str(s)]);
    tempBestPop = getFinalPop(bestPopListAll);
    [bestPop(:,s)] = tempBestPop;
    bestPopErrors(:,s) = bestPopValAll;
    
    for i = 1:length(tempBestPop)
        popListBest(1:i, popCnt) = sort(tempBestPop(1:i));
        popCnt = popCnt + 1;
    end
end
%%

stimType = 4;

currPop = bestPop(:,stimType);
popImprovement =  -diff(bestPopErrors(:, stimType))./(bestPopErrors(1:9, stimType))*100;
homoImprovement = [];
hetImprovement = [];
currPopInd = 2;
for currPopInd = 2:10
    if sum(bestPop(currPopInd, stimType) == bestPop(1:(currPopInd-1), stimType) > 0)
        homoImprovement(end+1) = popImprovement(currPopInd-1);
    else
        hetImprovement(end+1) = popImprovement(currPopInd-1);
    end
end
figure; plot(1, homoImprovement, 'r*', 2, hetImprovement, 'g*'); axis([0 3 0 max(popImprovement)+3])

ranksum(homoImprovement, hetImprovement)
%%
stimType = 1;

plot(1:10, bestPopErrors(:,stimType), 'k');
hold on;
currPopInd = 2;
for currPopInd = 2:10
    if sum(bestPop(currPopInd, stimType) == bestPop(1:(currPopInd-1), stimType) > 0)
        colorInd = 'r';
    else
        colorInd = 'g';
    end
    plot(currPopInd, bestPopErrors(currPopInd,stimType), ['*' colorInd]);
end
ylabel('Decoding error (rmse)');
xlabel('Population size');


for i = 1:2
    reconErrorPan(i).select();
    hold on;
    boundedline(1:numSearchIters , avgErr(:,i), avgStdErr(:,i), 'g', 'transparency', .2);
    plot(1:numSearchIters, bestPopValAll(i,:), '-b');
end
reconErrorPan.ylabel('Decoding error (rmse)');
reconErrorPan.xlabel('Population size');
reconErrorPan.fontsize = 6;

currPop = bestPop(:,stimType);
popImprovement =  -diff(bestPopErrors(:, stimType))./(bestPopErrors(1:9, stimType))*100;
homoImprovement = [];
hetImprovement = [];
currPopInd = 2;
for currPopInd = 2:10
    if sum(bestPop(currPopInd, stimType) == bestPop(1:(currPopInd-1), stimType) > 0)
        homoImprovement(end+1) = popImprovement(currPopInd-1);
    else
        hetImprovement(end+1) = popImprovement(currPopInd-1);
    end
end
figure; plot(1, homoImprovement, 'r*', 2, hetImprovement, 'g*'); axis([0 3 0 max(popImprovement)+3])

ranksum(homoImprovement, hetImprovement)