function [val] = meanAbsDiff(inMat)

% each column of inMat is a distinct vector
distance = 'minkowski';
minkVal = 1;

% first compute mean vector
meanVec = mean(inMat,2);
dists = zeros(1, size(inMat,2));

for i = 1:size(inMat,2)
    dists(i) = pdist([inMat(:,i) meanVec]', distance, minkVal);
end

val = mean(dists);