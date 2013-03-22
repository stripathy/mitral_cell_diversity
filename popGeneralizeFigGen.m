% population generalization figure - showing that hetero pops generalize to
% other stim better than homo pops

%% compute PSD's for example populations and example stim stats

% this is a file which contains 2 example pops (44 and 155) and 4 stim
% stats for 200 stimuli each (slen = 512)
cd('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\diffStimStats\finSims\sim10');
load('reconsMerged');
freqSimStimInds = [1 2 4 6 10];

[freq, stimSpect, reconSpect, infoSpect] = computeFreqSpectra(stimParms, reconStructOut, popAssignInds);
%recon and info spect are organized s.t rows are diff stim, and columns are
%diff populations
%% plot example stim and recons for a hetero and homo pop

% load data

cd('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\diffStimStats\finSims\finData\');
load reconPlusAnal2

%% plot example recons
colorCodePlot(1,:) = [1 0 0];
colorCodePlot(2,:) = [0 1 0];
colorCode = colormap(jet(length(keepCells)));

%% plot some example recons between 
currStimInd = 66;%20;
popInds = [44 155];
% popInds = [44 624];
for i = 1:size(popInds,2)
    plotInds(i) = find(currStimInd == popAssignInds.stimInd & popInds(i) == popAssignInds.popInd);
end

a = figure(); subplot(3,1,1:2); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 2) % axis([1 slen -3 3]);
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

currStimStat = 6;

stimInds = stimStatInds==currStimStat;
reconDiffVals = popValsAll(stimInds, popInds(1)) - popValsAll(stimInds, popInds(2));

[bestVal bestInds] = sort(reconDiffVals, 'descend');


%% set up global sturcture for figure
popInds = [44 155];
stimStatExamples = [1 2 6];
stimExampleInds = [7 16 45];
colorCode = colormap(jet(length(keepCells)));

%% plot the glm parms for a few pops
glmParmPan = panel();
glmParmPan.pack('v', [25 25 25 -1 ]);

for j = 1:length(popInds)
    currGlmPan = glmParmPan(j);
    currGlmPan.pack('h', [45 45 -1]);
    kFilt = currGlmPan(1);
    hFilt = currGlmPan(2);
    bFilt = currGlmPan(3);
    
currCellInds = pop(popInds(j)).dat;


% figure; 
% box off;
allCellInds = unique(nonzeros(currCellInds));
for i = 1:length(allCellInds)
    
    kFilt.select(); hold on; plot(-49:0, cellModels(allCellInds(i)).k, 'Color', colorCode(allCellInds(i),:) , 'LineWidth' , 1); 
    axis([-50 0 -.2 1.1]); xTickPos = [-50 -25 0]; hold on;
    set(gca, 'XTick', xTickPos);
    hFilt.select(); plot(cellModels(1).iht,log10(exp(cellModels(allCellInds(i)).ihbas*cellModels(allCellInds(i)).ih)), 'Color', colorCode(allCellInds(i),:), 'LineWidth' , 1); hold on;
    axis ([0 60 -1.6 2.2]);
    yTickPos = [-1:2]; yTickLabels = [.1 1 10 100];
    set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabels);
    bFilt.select(); plot(1, log10(exp(cellModels(allCellInds(i)).dc)),  '.', 'Color', colorCode(allCellInds(i),:)); hold on;
    axis ([0 2 -1.6 2.2]);
    set(gca, 'YTick', yTickPos, 'YTickLabel', [], 'XTick', []);
end
currGlmPan.margin = 5;
end
glmParmPan.margin = 5;
glmParmPan.fontsize = 6;



%plot all stim stats in frequency domain
stimFreqPan = glmParmPan(3);
stimFreqPan.pack('h', 2);
allFreqPan = stimFreqPan(1);

figure;
allFreqPan = panel();
plotStims = [1 2 6];

tapers=3; avg=1; fs = 1000;
pars = struct ('Fs', fs, ...
    'tapers', [tapers 2*tapers-1], ...
    'pad', 0, ...
    'trialave', 1);
for k = 1:numStimStats
    [S(:,k), f] = mtspectrumc(outSignal(:,1+(k-1)*numSignalsPerStat:k*numSignalsPerStat), pars);
end
allFreqPan.select();
plot(f,S(:,plotStims), 'LineWidth' , 1); ylabel('Power'); xlabel('Frequency (Hz)');
axis([0 100 0 .05]);
allFreqPan.fontsize = 6;
% subplot(1,5,1:2);box off;
% subplot(1,5,3:4); box off;
% subplot(1,5,5);box off;


%% plot stim and recons

allReconsPan = panel();
freqInds = 1:53;

allReconsPan.pack(length(stimStatExamples),1);
for k = 1:length(stimStatExamples)
    currStimStat = stimStatExamples(k);
    currStim = (currStimStat-1)*numSignalsPerStat + stimExampleInds(k);
        
    stimExample = allReconsPan(k,1);
    stimExample.pack('h', [1]);
    
    stimReconPan = stimExample(1);
%     freqPan = stimExample(2);
    
    stimReconPan.pack('v', [2/3 -1]);
    
    stimPan = stimReconPan(1);
    rasterPan = stimReconPan(2);
    
    
currStimInd = currStim;%20;
% popInds = [44 624];
for i = 1:size(popInds,2)
    plotInds(i) = find(currStimInd == popAssignInds.stimInd & popInds(i) == popAssignInds.popInd);
end

stimPan.select(); plot(1:slen, outSignal(:,currStimInd), 'k', 'LineWidth' , 1) % axis([1 slen -3 3]);
testSpikesPlot = [];
hold on;
rasterColors = [];
for i = 1:2%length(currCellInds)
    tempCurrCellInds = reconStructIn(plotInds(i)).popMakeup;
    plot(1:slen, reconStructOut(plotInds(i)).optStim, 'Color', colorCodePlot(i,:), 'LineWidth' , 1)
    testSpikesPlot = vectCat(testSpikesPlot, reconStructIn(plotInds(i)).testSpikes);
    
%     rasterColors(end+1:end+nnz(tempCurrCellInds),:) = repmat(colorCodePlot(i,:), nnz(tempCurrCellInds),1);
    rasterColors(end+1:end+nnz(tempCurrCellInds),:) = colorCode(tempCurrCellInds,:);
    if i ~= 2
        rasterColors(end+2,:) = [1,1,1]; %white
        testSpikesPlot = vectCat(testSpikesPlot, [0 0]);
    end
end
axis([0 slen -3 3]);
% legend('real stim', 'same recon1', 'same recon2','diff recon');
ylabel('Stimulus (a.u.)'); set(gca,'XTickLabel', []); axis tight; box off;
rasterPan.select(); createRaster3(testSpikesPlot',1, slen, rasterColors); %ylabel('Neuron ID');
xlabel('Time (ms)'); set(gca,'YTick', []);
set(gca,'TickDir', 'out')


%     plot freqPan;
%     currFreqInd = find(stimStatExamples(k) == freqSimStimInds);
%     freqPan.select();
%     plot(freq(freqInds), stimSpect(currFreqInd).dat(freqInds), 'k', 'LineWidth' , 1); hold on;
%     for i = 1:2%length(currCellInds)
%         plot(freq(freqInds), reconSpect(currFreqInd,i).dat(freqInds), 'Color', colorCodePlot(i,:), 'LineWidth' , 1);
%     end
%     axis tight;
    
    stimReconPan.de.margin = 2;
    stimExample.margin = 10;
end
allReconsPan.margin = 10;
allReconsPan.fontsize = 6;

%% plot stimuli population comparison panel
stimComps = [1 2; 1 6];
homoInds = 1:numModels;
heteroInds = numModels+1:numPops;

for i = 1:size(popAvgVal,2)
    for j = 1:numPops
        popPositionsStim(j,i) = find(sortedPopsStim(:,i) == j);
    end
end

figure;
compPan = panel();

compPan.pack(2, 2);
% stimCompExBothPan = compPan(1);

for i = 1:2
    stimInd1 = stimComps(i,1);
    stimInd2 = stimComps(i,2);
    
    compPan(1, i).select();
%     stimCompExPan(1,i).select();
    % compute average pairwise correlation of pop position across stim stats




% stimCompExPan.select(); 
% subplot(2,2,i);
hold on;
plot([1 numPops], [1 numPops], 'Color', [.4 .4 .4]);
plot(popPositionsStim(heteroInds,stimInd1), popPositionsStim(heteroInds,stimInd2),'g.')
plot(popPositionsStim(homoInds,stimInd1), popPositionsStim(homoInds,stimInd2),'r.')

plot(popPositionsStim(popInds(2),stimInd1), popPositionsStim(popInds(2),stimInd2),'g*')
plot(popPositionsStim(popInds(1),stimInd1), popPositionsStim(popInds(1),stimInd2),'r*')

set(gca,'YDir','reverse','XDir','reverse');
xlabel(['stim ', num2str(stimInd1)]);
ylabel(['stim ', num2str(stimInd2)]);
axis([0 244 0 244]);
set(gca, 'XTick', [1 100 200])
set(gca, 'YTick', [1 100 200])
end




% plot pairwise correlation panel for homo and hetero
% stimCorrPan = compPan(3);

numStimStats = stimParms.numStimStats;

bestPops = sortedPopsStim;

removeStims = [7 9];
includeStims = setdiff(1:numStimStats, removeStims);

homoCorr = 1 - pdist(popPositionsStim(homoInds,includeStims)', 'correlation');
heteroCorr = 1 - pdist(popPositionsStim(heteroInds,includeStims)', 'correlation');

homoSpread = pdist(popPositionsStim(homoInds,includeStims)', 'correlation');
heteroSpread = pdist(popPositionsStim(heteroInds,includeStims)', 'correlation');

homoCorrSq = squareform(homoCorr);
heteroCorrSq = squareform(heteroCorr);



% (1/(length(heteroInds)-1))*sum( (popPositionsStim(heteroInds, 1) - mean(1:numPops)).*(popPositionsStim(heteroInds, 6) - mean(1:numPops)))...
%     /(std(popPositionsStim(heteroInds, 1))*std(popPositionsStim(heteroInds, 6)));
% (1/(length(homoInds)-1))*sum( (popPositionsStim(homoInds, 1) - mean(1:numPops)).*(popPositionsStim(homoInds, 6) - mean(1:numPops)))...
%     /(std(popPositionsStim(homoInds, 1))*std(popPositionsStim(homoInds, 6)))

% heteroPts = popPositionsStim(heteroInds, 1);
% 
% dHet=abs(-popPositionsStim(heteroInds, 1)+popPositionsStim(heteroInds, 2))/sqrt(2);
% dHom=abs(-popPositionsStim(homoInds, 1)+popPositionsStim(heteroInds, 2))/sqrt(2);
% 
%     for i = 1:numel(heteroInds)
%         currPt = [popPositionsStim(heteroInds(i), 1); popPositionsStim(heteroInds(i), 2)];
%         d(i) = point_to_line(currPt, [0; 0], [1; 1]);
%     end
%     

% stimCorrPan.select();
compPan(2,2).select();
plot(homoSpread, heteroSpread, 'k.')
hold on;
plot([0 1], [0 1], 'Color', [.4 .4 .4])
xlabel('Homog- spread');
ylabel('Heterog- spread');

[pval, hyp] = signrank(heteroSpread, homoSpread);

compPan.de.margin = 10;
compPan.fontsize = 6;

%% plot average position for het and hom across stims
meanPosHom = mean(mean(popPositionsStim(homoInds,includeStims)));
stdErrPosHom = myStdErr(nonzeros(popPositionsStim(homoInds,includeStims)));

meanPosHet = mean(mean(popPositionsStim(heteroInds,includeStims)));
stdErrPosHet = myStdErr(nonzeros(popPositionsStim(heteroInds,includeStims)));

compPan(2,1).select();
hold on;
errorbar([meanPosHom meanPosHet], [stdErrPosHom stdErrPosHet], 'Color', 'k', 'LineStyle', 'none');
bar([meanPosHom meanPosHet], 'baseValue', numPops);
set(gca,'YDir','reverse')
ylabel('Mean population position');
currxTickLabel = {'Hom-', 'Het-'};
set(gca,'XTick', [1 2], 'XTickLabel', currxTickLabel)
axis([0 3 50 numPops]);

[pval, hyp] = signrank(mean(popPositionsStim(homoInds,includeStims)), mean(popPositionsStim(heteroInds,includeStims)));
% 
% plot(ones(10,1), mean(popPositionsStim(homoInds,:)), 'r.'); hold on;
% plot(ones(10,1)+1, mean(popPositionsStim(heteroInds,:)), 'g.')
% 
% plot([ones(10,1) ones(10,1)+1]', [mean(popPositionsStim(homoInds,:))', mean(popPositionsStim(heteroInds,:))']', 'k')
%% plot frequency panel
p = panel();

freqInds = 1:53;
p.pack(length(stimStatExamples),1);
for k = 1:length(stimStatExamples)
    currFreqInd = find(stimStatExamples(k) == freqSimStimInds);
    freqPan = p(k,1);
    freqPan.select();
    plot(freq(freqInds), stimSpect(currFreqInd).dat(freqInds), 'k'); hold on;
    for i = 1:2%length(currCellInds)
        plot(freq(freqInds), reconSpect(currFreqInd,i).dat(freqInds), 'Color', colorCodePlot(i,:));
    end
    axis tight;
end

%% make entire figure
close all
clf
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0 0 18.3 12.7]);

popGenFig = panel();
popGenFig.pack('h', [35 40 25]);

glmParmPan = popGenFig(1);
allReconsPan = popGenFig(2);
compPan = popGenFig(3);




    
