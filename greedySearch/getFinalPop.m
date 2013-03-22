function [addedNeuron] = getFinalPop(bestPopListAll)
% figure out neuron addition path
currBestPop = [];
addedNeuron = [];

numSearchIters = length(bestPopListAll);
for i = 1:numSearchIters
    addedNeuron(i) = multisetdiff(bestPopListAll(i).dat, currBestPop);
    currBestPop = bestPopListAll(i).dat;
end
finBestPop = currBestPop;
for i = 1:numSearchIters
    numNeuronsInBest(i) = sum(addedNeuron(i)==finBestPop);
end