function [optStim, exitflag, fval, hessian] = bayesStimDecoder(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat)
%inputs: cellModels is a vector of structs where each element of the struct
%is as passed back from the Pillow LNP fit algorithm
% spikeTimes is a matrix where each column is the relevant decoding spike
% times for the ith neuron

%bin input spike times
tPts = stimInterval(1):dtStim:stimInterval(2);
spikeMat = histc(spikeTimes, tPts);
slen = length(tPts);

numCells = length(cellModels);

%initialize matrices for kfilts and nonkfilts
staWin = length(cellModels(1).k);
kFiltMat = zeros(staWin, numCells);
hFiltFact = zeros(slen, numCells);

for i = 1:numCells
    kFiltMat(:,i) = cellModels(i).k;
    
    ihTemp = cellModels(i).ihbas*cellModels(i).ih;
%     hFiltTemp = [0; downsample(ihTemp, 1/cellModels(i).dt)];
    hFiltTemp = [downsample(ihTemp, 1/cellModels(i).dt)];
    tempConv = conv(spikeMat(:,i), hFiltTemp) + cellModels(i).dc;
    hFiltFact(1:slen,i) = tempConv(1:slen);
end

priorCovMatDet = det(stimCovMat);
priorInvCovMat = inv(stimCovMat);

% i should write something that computes the first estimated stim via opt
% lin decoding

% currEstStim = randn(slen,1)*stimCovMat(1,1); %this makes a random GWN stim sample
currEstStim = getNoisyStim(slen, 1, 3);
% [tsp,Vmem,Ispk] = simGLM(cellModels(i),currEstStim);

inputParms.spikeMat = spikeMat;
inputParms.dtStim = dtStim;
inputParms.kFiltMat = kFiltMat;
inputParms.hFiltFact = hFiltFact;

inputParms.priorCovMatDet = priorCovMatDet;
inputParms.priorInvCovMat = priorInvCovMat;

% opts = optimset('Gradobj','on','Hessian','on', 'TolFun', 1e-12, 'TolX', 1e-12, 'maxIter', 1000, 'Display', 'on');
opts = optimset('Gradobj','on','Hessian','on', 'maxIter', 1000, 'TolFun', 1e-8, 'TolX', 1e-8, 'Display', 'off');

prs0 = currEstStim;
fxn = @(prs)bayesDecoder_logli(prs,inputParms);

exitflag = 0;
maxNumRuns = 1;
numRuns = 0;
while exitflag == 0 && numRuns < maxNumRuns
    [prs,fval,exitflag] = fminunc(fxn,prs0,opts);
    prs0 = prs;
    numRuns = numRuns + 1;
end

optStim = prs;
if nargout > 2 % Compute Hessian if desired
    [fval,gradval,hessian] = bayesDecoder_logli(prs, inputParms);
end
% figure; plot(1:slen, currEstStim, 1:slen, optStim);