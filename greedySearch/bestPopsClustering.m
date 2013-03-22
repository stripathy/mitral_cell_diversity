% try to do some clustering of the best pops

%represent best populations via a 1-44 neuron count representation
popSizeCon = 5; % largest pop size to consider

bestPopCntRep = zeros(44, size(bestPop,2));
for i = 1:size(bestPop,2)
    bestPopCntRep(:,i) = histc(bestPop(1:popSizeCon,i), 1:44);
end
    
numRandomPops = 25;
randPopCntRep = zeros(44, numRandomPops);

clear randPopList
numBestPops = size(popListBest,2);
parfor i = 1:numRandomPops
%     popSize(i) = randi(numSearchIters);
    randPopList(:,i) = randsample(1:44, popSizeCon, 1);    
    randPopCntRep(:,i) = histc(randPopList(:,i), 1:44);
    
%     [stimDiffRand(:,i) histDiffRand(:,i) biasDiffRand(:,i)] = glmDiff(cellModels(tempPopList));
end

dists = pdist([bestPopCntRep(:,1:50) randPopCntRep]');
figure; imagesc(squareform(dists))

clusterdata([bestPopCntRep randPopCntRep], 2);

figure;
Z = linkage([bestPopCntRep(:,1:50) randPopCntRep]', 'average');
dendrogram(Z, 0, 'labels', dendLabs)

clear dendLabs
cnt = 1;
for i = 1:numStimStats
    for j = 1:5
        dendLabs{cnt} = num2str(i);
        cnt = cnt + 1;
    end
end
for i = 1:numRandomPops
    dendLabs{cnt} = num2str(0);
    cnt = cnt + 1;
end

eucDistMatShape = squareform(dists);

options.dims = 1:2;
[Y, R, E] = Isomap(eucDistMatShape, 'k', 100, options); 
figure
hold on;
for i = 1:75%length(allRatesVec)
plot(Y.coords{1}(1,i),Y.coords{2}(2,i), '.')
text(Y.coords{1}(1,i)+.05,Y.coords{2}(2,i), num2str(i));
end