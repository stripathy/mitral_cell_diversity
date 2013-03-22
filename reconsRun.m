function [] = reconsRun(splitInds)
% script for running lots of decoding

% matlabpool local;

% cd('/home/stripathy/matlab/stim_recon/diffCorrs/sim3');
cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\naturalComp\sim2');

for i = 1:length(splitInds)
    
    load(['split', num2str(splitInds(i))]);
    
%     [reconStructOut, decodingTime] = reconRunfile(cellModels,currReconStructIn, stimCovMat, outSignal);
    [reconStructOut, decodingTime] = reconRunfile(cellModels,currReconStructIn, stimCovMat, outSignal);
    
    save (['split', num2str(splitInds(i)), 'Out']);
    clear reconStructOut
    
end

matlabpool close;
