function [] = reconsRunTemp(splitInds)
% script for running lots of decoding

matlabpool local;

cd('/home/stripathy/matlab/stim_recon/diffCorrs/sim4');
% cd ('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\diffStimStats\sim8');

for i = 1:length(splitInds)
    
    load(['split', num2str(splitInds(i))]);
    
%     [reconStructOut, decodingTime] = reconRunfile(cellModels,currReconStructIn, stimCovMat, outSignal);
    [reconStructOut, decodingTime] = reconRunfile(cellModels,currReconStructIn, stimCovMat, outSignal);
    
    save (['split', num2str(splitInds(i)), 'Out']);
    clear reconStructOut
    
end

matlabpool close;
