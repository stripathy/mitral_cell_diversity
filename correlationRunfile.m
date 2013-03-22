

%% generate all pairwise crosscorrelations

% This code assumes that there is a single long stimulus and multiple (even
% number) of repeats

%get all spikes from all neurons
allSpikes = []; trialsPerNeuron = [];
for i = 1:2%length(keepCells)
    allSpikes = vectCat(allSpikes, simSpikes(i).dat);
end
trialsPerNeuron(1:i) = numSimRepeats;

% run xcorr algorithm
tic
[xcorrs, acorrs1, acorrs2, xcorrsPairwise, autoCorrs]=crossCorrAllCells(allSpikes', trialsPerNeuron, slen);
toc