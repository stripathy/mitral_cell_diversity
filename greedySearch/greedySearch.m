function [] = greedySearch(objFxnInd, stimInd)

%% code to run greedy search on best pops vs stims
tic
% load data - stim and spike trains over neurons
load(['template', num2str(stimInd)]);
% if stimInd == 1
%     load templateHigh
% elseif stimInd == 2
%     load templateLow
% elseif stimInd == 3
%     load templateHighLow
% else
%     load templateSine
% end
numSearchIters = 10;
bestPopThresh = 1;

[simSpikes, simSpikesInds] = genSpikesFromStim(outSignal, stimParms, cellModels);

% generate a pop list - this stores all populations simulated
pop = struct;
popOld = struct;

% if objFxnInd == 2
%     for i = 1:size(stimCovMat,3)
%         stimEnt(i) = .5*sum(myLog2(eig(stimCovMat))) + .5*slen*myLog2(2*pi*exp(1));
%     end
% end

for m = 1:numSearchIters
    
    % make populations to search
    [popNew, newPopList] = greedySearchMakePops(popOld, neuronsToSim, m);
    numPops = length(popNew);
    
    % code to assign spikes to pops and make reconStructIn
    
    [reconStructIn, popAssignInds] = assignPopSpikes(popNew, simSpikes, simSpikesInds, stimParms);
    
    % decode stims from reconStructIn
    % but reconStructIn traditionally gets broken up first
    [reconStructOut] = reconRunfile(cellModels,reconStructIn, stimCovMat, outSignal);
    
    % find the best pops given the objective function (rmse, max info/sec)
    if objFxnInd == 1 || objFxnInd == 4
        [reconErrorPopAvg, reconErrorPopAll] = objFxnGreedySearch(reconStructOut, popAssignInds, stimParms, objFxnInd, outSignal);
        if objFxnInd == 1
            [bestVals bestPopsAllStims] = sort(mean(reconErrorPopAvg,1), 'ascend');
        else %objFxnInd = 4 - take best ranks acros stim
            [sortedVals, popRanks] = sort(reconErrorPopAvg, 2, 'ascend');
            [bestVals bestPopsAllStims] = sort(mean(popRanks,1), 'ascend');
        end
            
    elseif objFxnInd == 2 || objFxnInd == 3
        [respEntPopAvg] = objFxnGreedySearch(reconStructOut, popAssignInds, stimParms, objFxnInd);
        totalEnt = repmat(stimEnt', 1, numPops) - (respEntPopAvg  + .5*slen*myLog2(2*pi*exp(1)));
        totalEntPerSec = totalEnt/(slen*.001);
        [bestVals bestPopsAllStims] = sort(mean(totalEntPerSec,1), 'descend');
        
        if objFxnInd == 3
            for i = 1:numStimStats
                for j = 1:numPops
                    popSumFR(i,j) = sum(meanFRPerStim(newPopList(:,j),i));
                end
            end
            totalEntPerSpike = totalEntPerSec./popSumFR;
            [bestVals bestPopsAllStims] = sort(mean(totalEntPerSpike,1), 'descend');
        end
            
    end
    
%     [sorted bestPopsAllStims] = sort(mean(popAvgVal));  
    bestPopId = bestPopsAllStims(1:bestPopThresh);
    bestPopList = newPopList(:,bestPopsAllStims(1:bestPopThresh))
    bestPopVal = bestVals(1:bestPopThresh);
    
%     
%     for i = 1:numStimStats+1
%         bestPopId(:,i) = bestPopsAllStims(i,1:bestPopThresh);
%         bestPopList(:,:,i) = newPopList(:,bestPopsAllStims(i,1:bestPopThresh));
%     end
    
    bestPopIdAll(m).dat = bestPopId;
    bestPopListAll(m).dat = bestPopList;
    bestPopValAll(m) = bestPopVal;
    popValsAll(m).dat = reconErrorPopAll;
%     reconsAll(m).in = reconStructIn;
%     reconsAll(m).out = reconStructOut;
    
    popOld = popNew(unique(bestPopId));
    % save some data
%     save('greedySearchIter' [m]);
    
    %delete some data
    clear reconStructIn reconStructOut bestPopId bestPopList bestPopVal
        
    %run next iteration of greedy search;
    
end

save (['greedySearch', num2str(objFxnInd), num2str(stimInd)]);

t = toc