% best pops why plots

%% illustrate the problem of asking what makes certain populations good

% popInds =  [822 624 424 838];
popInds =  [155];
numRows = length(popInds);
figure;
for j = 1:length(popInds)
    currCellInds = popList(:,popInds(j));
    % currCellInds = [10 27 33 41 44];
    % currCellInds = [27 33];
    colorCode = colormap(jet(length(keepCells)));
    
    
    allCellInds = unique(nonzeros(currCellInds));
    for i = 1:length(allCellInds)
        
        subplot(numRows,5,(1:2) + (j-1)*5);plot(-49:0, cellModels(allCellInds(i)).k, 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 2); axis([-50 0 -.5 1.2]); hold on;
        ylabel('Stim units (a.u.)'); xlabel('Time (ms)');
        subplot(numRows,5,(3:4)+ (j-1)*5);plot(cellModels(1).iht,log10(exp(cellModels(allCellInds(i)).ihbas*cellModels(allCellInds(i)).ih)), 'Color', colorCode(allCellInds(i),:), 'LineWidth' , 2); hold on;
        axis ([0 60 -1.6 2.2]);
        yTickPos = [-1:2]; yTickLabels = [.1 1 10 100];
        set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabels);
        xlabel('Time (ms)');
        ylabel('Gain');
        subplot(numRows,5,5+ (j-1)*5);plot(1, log10(exp(cellModels(allCellInds(i)).dc)),  '*', 'Color', colorCode(allCellInds(i),:)); hold on;
        axis ([0 2 -1.6 2.2]);
        set(gca, 'YTick', yTickPos, 'YTickLabel', [], 'XTick', []);
    end
end

%% illustrate regression approach to problem

% see how stim 1&3 depends on predictor variables

% comparison to FR alone
stimInd = 3;
for i = 1:numPops
    popMeanFR(i) = mean(meanFRPerStim(popList(:,i),stimInd));
end

X = [ones(numPops, 1) popMeanFR'];
y = reconRMSEMean(stimInd,:)';
betahat = inv(X'*X)*X'*y;
resids = y - X*betahat;

figure;
plot(popMeanFR, reconRMSEMean(stimInd,:),'.' , popMeanFR, X*betahat, 'k.')
xlabel('mean FR (Hz)'); ylabel('reconstruction error (rmse)');

figure;
plot(popMeanFR, resids,'.');
xlabel('mean FR (Hz)'); ylabel('residuals');

% comparison to FR and meanStimDiff alone
stimInd = 1;
depVar = reconRMSEMean(stimInd,:)';
for i = 1:numPops
    popMeanFR(i) = mean(meanFRPerStim(popList(:,i),stimInd));
end

indepVar1 = popMeanFR;
indepVar2 = meanStimDiff;

% % color cells directly by FR
colorCodeTotal = colormap(hot(max(floor(depVar*1000)) - min(floor(depVar*1000)) + 1));
for i = 1:numPops
    colorCode(i,:) = colorCodeTotal(floor(depVar(i)*1000) - min(floor(depVar*1000)) +1,:);
end

figure;
hold on;
for i = 1:numPops
    plot(indepVar1(i), indepVar2(i), '.', 'Color', colorCode(i,:));
end
colormap hot; caxis([min(depVar) max(depVar)]); colorbar;
xlabel('mean FR (Hz)'); ylabel('mean stim diff (L1)');
axis tight;




