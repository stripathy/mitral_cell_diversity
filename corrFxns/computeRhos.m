function [rhos] = computeRhos(xcorrs, acorrs1, acorrs2)

windowLen = 5;
corrFxnLen = size(xcorrs,2);
windowFun = (windowLen - abs(-(windowLen-1)/2:(windowLen-1)/2))/windowLen;
windowFun = windowFun';

rhos = zeros(size(xcorrs,1),1);
corrInds = (corrFxnLen-1)/2 + 1 +((-(windowLen-1)/2):((windowLen-1)/2)); 
for i = 1:size(xcorrs,1)
    rhos(i) = xcorrs(i,corrInds)*windowFun/sqrt(acorrs1(i,corrInds)*windowFun*acorrs2(i,corrInds)*windowFun);
end