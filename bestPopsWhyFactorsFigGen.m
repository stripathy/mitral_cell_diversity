% best pops why plots

%% 
clear projFilters
projFilters = [pMatrixA(:,1:5)];
% 
% 
% projFilters = [pMatrixK(:,1:numEigsK) pMatrixH(:,1:numEigsH) pMatrixB(:,1:numEigsB)];
% 
% projFilterMeans = mean(projFilters);
% projFilterStds = std(projFilters);
% 
% projFilters = zscore(projFilters);
% 
% projFilters = [pMatrixA(:,1:5)];


%% compute charateristics of each pop along each dimenson

clear *FiltParms
for i = 1:numPops
    for j = 1:maxPopSize
        filtParmInds(j) = find(popList(j,i) == neuronsToSim);
    end
    meanProjFiltParms(i,:) = mean(projFilters(filtParmInds,:),1);
    stdProjFiltParms(i,:) = std(projFilters(filtParmInds,:),1);
    if numel(unique(filtParmInds)) > 1
        skewProjFiltParms(i,:) = skewness(projFilters(filtParmInds,:),[],1);
        kurtProjFiltParms(i,:) = kurtosis(projFilters(filtParmInds,:),[],1);
    else
        skewProjFiltParms(i,:) = zeros(1, size(projFilters,2));
        kurtProjFiltParms(i,:) = zeros(1, size(projFilters,2));
    end
end

%% define popAvgVals as recon accuracy, not error
popAvgVal = 1 - popAvgVal;
%% plot a couple examples illustrating regression approach, i.e. rmse of stim stat 1 vs glm parm pca1

validPopInds = 45:244;
 
% regressExFig = panel();
% regressExFig.select();
% regressExFig.fontsize = 6;
stimStats = [1 2];
figure;
for i = 1:length(stimStats)
    
subplot(2,2,i); plot(meanProjFiltParms(validPopInds,1), popAvgVal(validPopInds,stimStats(i)),'.'); 
axis tight; box off;
set(gca, 'FontSize', 7);
subplot(2,2,i+2); plot(stdProjFiltParms(validPopInds,3), popAvgVal(validPopInds,stimStats(i)),'.')
axis tight; box off;
set(gca, 'FontSize', 7);
end
% regressExFig.fontsize = 6;

% plot(stdProjFiltParms(validPopInds,3), popAvgVal(validPopInds,1),'.')

%% plot a couple examples illustrating regression approach, i.e. rmse of stim stat 1 vs glm parm pca1

validPopInds = 45:244;
 
% regressExFig = panel();
% regressExFig.select();
% regressExFig.fontsize = 6;
stimStats = [1];
figure;
for i = 1:length(stimStats)
    yDat = popAvgVal(validPopInds,stimStats(i));
    xDat = meanProjFiltParms(validPopInds,1);
    subplot(2,1,1); plot(xDat, yDat, '.');
    hold on;
    linFit = polyfit(xDat, yDat, 1);
    yFit = polyval(linFit, xDat);
    plot(xDat, yFit, 'k');
    axis tight; box off;
    set(gca, 'FontSize', 7);
    rsqEx(1) = corr(yDat, yFit).^2;
    
    xDat = stdProjFiltParms(validPopInds,3);
    subplot(2,1,2); plot(xDat, yDat, '.');
    hold on;
    linFit = polyfit(xDat, yDat, 1);
    yFit = polyval(linFit, xDat);
    plot(xDat, yFit, 'k');
    axis tight; box off;
    set(gca, 'FontSize', 7);
end
% regressExFig.fontsize = 6;

% plot(stdProjFiltParms(validPopInds,3), popAvgVal(validPopInds,1),'.')

%% compute beta coeffs and R^2 values, on pooled GLM PCA space

% compute R^2 values for different kinds and amounts of variables

% predictors (X's) to test: 
% 1. which PCAs to use? 1-3
% 2. how many PCA population moments to use (1-3)? population means, stds,
% skews?
figure;
for k = 1:length(stimStats)
    currStimStat = stimStats(k);
    
    numPCAs = 5;
    popMoments = 4;
    clear allPopParms rsqVal betaCoeffs
    allPopParms(:,:,1) = meanProjFiltParms;
    allPopParms(:,:,2) = stdProjFiltParms;
    allPopParms(:,:,3) = skewProjFiltParms;
    allPopParms(:,:,4) = kurtProjFiltParms;
    
    % validPopInds = 45:244;
    Y = zscore(popAvgVal(validPopInds,currStimStat));
    for i = 1:numPCAs
        for j = 1:popMoments
            regressors = [allPopParms(validPopInds,1:i, 1:j)];
            regressors = reshape(regressors, numel(validPopInds), i*j);
            
            X = [ones(numel(validPopInds), 1) zscore(regressors)];
            [b,bint,r,rint,stats] = regress(Y,X);
            
            
            rsqVal(i,j) = stats(1);
            betaCoeffs(i,j).dat = b;
            betaInt(i,j).dat = bint;
            
        end
    end
    rsqValKeep(k) = rsqVal(3,2);
    
    usePCAs = 3;
    
    subplot(1,length(stimStats),k);
    imagesc(rsqVal)
    caxis([0 1]);
%     colorbar();
    colormap('hot');
    set(gca, 'FontSize', 7);
    xlabel('number population moments');
    ylabel('number PCAs used');
%     zlabel('R^2 value');
end

k = 6; n = numel(validPopInds);
rsq = rsqValKeep(1);
rsqErrors(1) = sqrt((4*rsq*(1-rsq)^2*(n-k-1)^2)/((n^2-1)*(n+3)));
rsq = rsqValKeep(2);
rsqErrors(2) = sqrt((4*rsq*(1-rsq)^2*(n-k-1)^2)/((n^2-1)*(n+3)));

figure;
hold on;
bar([1, 2], rsqValKeep);
errorbar([1, 2], rsqValKeep, rsqErrors, 'k.')
ylabel('Percent variance explained (R^2)');




figure;
testStims = [1 2 8 4];

clear betaCoeffs betaInt
for i = 1:length(testStims)
    subplot(1, 2, i);
    Y = zscore(popAvgVal(validPopInds,testStims(i)));
%     
%             regressors = [allPopParms(validPopInds,1:i, 1:j)];
%         regressors = reshape(regressors, numel(validPopInds), i*j);
% %         
%         X = [ones(numel(validPopInds), 1) zscore(regressors)];
%         [b,bint,r,rint,stats] = regress(Y,X);
        
    X = [ones(numel(validPopInds), 1) zscore(meanProjFiltParms(validPopInds,1:usePCAs)) zscore(stdProjFiltParms(validPopInds,1:usePCAs))];
    [b,bint,r,rint,stats] = regress(Y,X);
    
    bar(b(2:end), 'g');
    betaCoeffs(i).dat = b;
    hold on; box off;
    errorbar(1:length(b(2:end)), b(2:end), b(2:end) - bint(2:end,1), b(2:end) - bint(2:end,2), '.k')
    betaInt(i).dat = bint;
    axis([0 size(X,2) -.5 1]);
    set(gca, 'FontSize', 7);
    set(gca, 'XTick', []);
    
    ylabel('Standardized regression coeffs');
end

%% plot beta coeffs for particular stimuli only
testStims = [1:10];

clear betaCoeffs betaInt
for i = 1:length(testStims)

    Y = zscore(popAvgVal(validPopInds,testStims(i)));
%     
%             regressors = [allPopParms(validPopInds,1:i, 1:j)];
%         regressors = reshape(regressors, numel(validPopInds), i*j);
% %         
%         X = [ones(numel(validPopInds), 1) zscore(regressors)];
%         [b,bint,r,rint,stats] = regress(Y,X);
        
    X = [ones(numel(validPopInds), 1) zscore(meanProjFiltParms(validPopInds,1:usePCAs)) zscore(stdProjFiltParms(validPopInds,1:usePCAs))];
    [b,bint,r,rint,stats] = regress(Y,X);
    betaCoeffs(i).dat = b;
    betaInt(i).dat = bint;
    rsqVals(i) = stats(1);

end

for i = 1:10
    betaCoeffsAll(:,i) = betaCoeffs(i).dat;
end

figure;

freqOrd = [2 8 1 4]; % for these stimuli
colormap redblue;
imagesc(betaCoeffsAll(2:end,freqOrd)); 
caxis([-1 1]);
    set(gca, 'FontSize', 7);
    colorbar;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'CTick', [-1:.5:1]);
    hcb = colorbar('YTickLabel',[-1:.5:1]);
set(hcb,'YTickMode','manual')

% 6 5 3 10];
