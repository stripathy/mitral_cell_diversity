[popAvgVal, popValsAll] = computeAvgErrors(reconStructOut, popAssignInds, stimParms, outSignal);

[blah, sortedPops] = sort(a);

[blah, sortedPopsStim] = sort(popAvgVal);

% compute average pairwise correlation of pop position across stim stats
for i = 1:size(popAvgVal,2)
    for j = 1:numPops
        popPositionsStim(j,i) = find(sortedPopsStim(:,i) == j);
    end
end

homoInds = 1:numModels;
heteroInds = numModels+1:numPops;

stimInd1 = 1;
stimInd2 = 10;

figure; 
hold on; plot(popPositionsStim(heteroInds,stimInd1), popPositionsStim(heteroInds,stimInd2),'g.')
plot(popPositionsStim(homoInds,stimInd1), popPositionsStim(homoInds,stimInd2),'r.')

set(gca,'YDir','reverse','XDir','reverse');
xlabel(['stim ', num2str(stimInd1)]);
ylabel(['stim ', num2str(stimInd2)]);


bestPops = sortedPopsStim;

sortedPopPositions = zeros(size(bestPops));
for i = 1:numStimStats
    for j = 1:numPops
        sortedPopPositions(bestPops(i,j), i) = j;
    end
end

removeStims = [];
includeStims = setdiff(1:numStimStats, removeStims);

homoCorr = 1 - pdist(popPositionsStim(homoInds,includeStims)', 'correlation');
heteroCorr = 1 - pdist(popPositionsStim(heteroInds,includeStims)', 'correlation');

homoCorrSq = squareform(homoCorr);
heteroCorrSq = squareform(heteroCorr);



% homoCorrIncSq = homoCorrSq(includeStims,includeStims);
% heteroCorrIncSq = heteroCorrSq(includeStims,includeStims);

figure;
plot(homoCorr, heteroCorr, '.', [0 1], [0 1])
% hold on;
% plot(mean(homoCorr), mean(heteroCorr), '*');

%% plot example recons
colorCodePlot(1,:) = [1 0 0];
colorCodePlot(2,:) = [0 1 0];
colorCode = colormap(jet(length(keepCells)));

%% plot some example recons between 
currStimInd = 7;%20;
popInds = [44 155];
% popInds = [44 624];
for i = 1:size(popInds,2)
    plotInds(i) = find(currStimInd == popAssignInds.stimInd & popInds(i) == popAssignInds.popInd);
end

figure; subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2) % axis([1 slen -3 3]);
testSpikesPlot = [];
hold on;
rasterColors = [];
for i = 1:2%length(currCellInds)
    tempCurrCellInds = reconStructIn(plotInds(i)).popMakeup;
    subplot(3,1,1:2); plot(1:slen, reconStructOut(plotInds(i)).optStim, 'Color', colorCodePlot(i,:), 'LineWidth' , 2)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    
%     rasterColors(end+1:end+nnz(tempCurrCellInds),:) = repmat(colorCodePlot(i,:), nnz(tempCurrCellInds),1);
    rasterColors(end+1:end+nnz(tempCurrCellInds),:) = colorCode(tempCurrCellInds,:);
    if i ~= 2
        rasterColors(end+1,:) = [1,1,1]; %white
        testSpikesPlot = vectCat(testSpikesPlot, 0);
    end
end
axis([0 slen -3 3]);
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('Stimulus (a.u.)'); set(gca,'XTickLabel', []); axis tight; box off;
subplot(3,1,3);createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('Time (ms)'); set(gca,'YTick', []);

%% find stim where diff between pop 1 and 2 is largest
cnt = 1;
for i = 1:stimParms.numStimStats
    for j = 1:numSignalsPerStat
        stimStatInds(cnt) = i;
        cnt = cnt + 1;
    end
end

currStimStat = 1;

stimInds = stimStatInds==currStimStat;
reconDiffVals = popValsAll(stimInds, popInds(1)) - popValsAll(stimInds, popInds(2));

[bestVal bestInds] = sort(reconDiffVals, 'descend');

