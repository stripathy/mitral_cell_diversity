% this is the version of bayesStimDecoder which accepts neurons which are
% coupled and described via the ggsim structure in Pillow's code
% this assumes each neuron in the pop recieves the same stim

function [optStim, exitflag, fval, hessian] = bayesStimDecoderCpl(ggsim, spikeTimes, stimInterval, dtDecoding, stimCovMat, initStim)

numCells = ggsim.numCells;
tPts = dtDecoding:dtDecoding:length(initStim)*dtDecoding;
% tPts = stimInterval(1):dtDecoding:stimInterval(2);
slen = length(tPts);
kMat = zeros(size(ggsim.k,1)/dtDecoding, numCells); %columns are the kfilter for each cell

% get kMat and spikeVecMat

spikeVecMat = zeros(slen, numCells);
spikeIndsMat = zeros(nnz(spikeTimes),numCells);
for i = 1:numCells
    spikeInds = ceil((nonzeros(spikeTimes(:,i))-dtDecoding*.001)/dtDecoding);
    spikeVec = zeros(slen, 1);
    spikeVec(spikeInds) = 1;
    
    %remember that kFilts end at t = 0
    
    kMat(:,i) = fastinterp2(ggsim.k(:,1,i),round(1/dtDecoding)) * dtDecoding; %this will need to get rescaled to match dtDecoding
    spikeVecMat(:,i) = spikeVec;
    spikeIndsMat(1:length(spikeInds),i) = spikeInds;
end

%compute histTermMat
histTermMat = zeros(length(tPts), numCells);
for i = 1:numCells
    hTermKernelMat = squeeze(ggsim.ih(:,i,:));
    hTermKernelMat = downsample(hTermKernelMat,round(dtDecoding/ggsim.dt));
    hTermKernelMat(1,:) = 0;
    
    for j = 1:numCells
        histTermMatTemp = spikeconv_mex(nonzeros(spikeIndsMat(:,j)),hTermKernelMat(:,j),[1,slen]);
        histTermMat(:,i) = histTermMat(:,i) + histTermMatTemp(1:slen); % this is the unchanging hist comp + dc comp
    end
    histTermMat(:,i) = histTermMat(:,i) + ggsim.dc(i);
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
opts = optimset('Gradobj','on','Hessian','on', 'MaxFunEvals', 5000000, 'MaxIter', 5000);
[prs,fval,exitflag] = fminunc(fxn,prs0,opts);

if nargout > 3 % Compute Hessian if desired
    [fval,gradval,hessian] = bayesStimDecoderLogli(prs, inputParms);
end

optStim = prs;
