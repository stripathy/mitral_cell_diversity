% this is the basic bayesStimDecoder, it seems to perform decoding
% correctly as long as the dt of the stimulus is 1 msec

function [optStim, exitflag, fval, hessian] = bayesStimDecoder1(cellModels, spikeTimes, stimInterval, dtDecoding, stimCovMat, initStim)

numCells = length(cellModels);
tPts = dtDecoding:dtDecoding:length(initStim)*dtDecoding;
% tPts = stimInterval(1):dtDecoding:stimInterval(2);
slen = length(tPts);
kMat = zeros(length(cellModels(1).k)/dtDecoding, numCells); %columns are the kfilter for each cell

histTermMat = zeros(length(tPts), numCells);
spikeVecMat = zeros(slen, numCells);
for i = 1:numCells
    spikeInds = ceil((nonzeros(spikeTimes(:,i))-dtDecoding*.001)/dtDecoding);
    spikeVec = zeros(slen, 1);
    spikeVec(spikeInds) = 1;
    
    %remember that kFilts end at t = 0
    
    kMat(:,i) = fastinterp2(cellModels(i).k,round(1/dtDecoding)) * dtDecoding; %this will need to get rescaled to match dtDecoding
    hTermKernel = cellModels(i).ihbas*cellModels(i).ih;
%     hTermKernel = cellModels(i).ih;
    hTermKernel = downsample(hTermKernel,round(dtDecoding/cellModels(i).dt));
    hTermKernel(1) = 0;
    
    histTermMatTemp = spikeconv_mex(spikeInds,hTermKernel,[1,slen]);
    histTermMat(:,i) = histTermMatTemp(1:slen) + cellModels(i).dc; % this is the unchanging hist comp + dc comp
    
    spikeVecMat(:,i) = spikeVec;
end

% currEstStim = randn(slen,1);
% stimCovMat = eye(slen);

inputParms.kMat = kMat;
inputParms.histTermMat = histTermMat;
inputParms.stimCovMatInv = inv(stimCovMat);
inputParms.stimCovMatDet = det(stimCovMat);
inputParms.spikeVecMat = spikeVecMat;
inputParms.dt = dtDecoding;

if nargin > 5
    currEstStim = initStim;
else
    currEstStim = getNoisyStim(slen*dtDecoding,dtDecoding, 3);
end

prs0 = currEstStim;
fxn = @(prs)bayesStimDecoderLogli(prs,inputParms);

% opts = optimset('Gradobj','on','Hessian','on', 'maxIter', 1000, 'TolFun', 1e-8, 'TolX', 1e-8);
% opts = optimset('Gradobj','on','Hessian','on', 'display', 'iter');
opts = optimset('Gradobj','on','Hessian','on', 'MaxFunEvals', 5000000, 'MaxIter', 5000, 'Display', 'off');
[prs,fval,exitflag] = fminunc(fxn,prs0,opts);

if nargout > 3 % Compute Hessian if desired
    [fval,gradval,hessian] = bayesStimDecoderLogli(prs, inputParms);
end

optStim = prs; 
