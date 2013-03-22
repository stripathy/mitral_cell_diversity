function [logli, grad, hessian] = logli_GLM(cellModels, spikeTimes, stim, stimCovMatInv, dtDecoding)
%computes the positive likelihood for a GLM model(s) given spikeTrains from
%each neuron and a stimulus. gradients and hessians are computed with
%respect to the stimulus vector

if nargin < 5
    dtDecoding = 1;
end

numCells = length(cellModels);
tPts = dtDecoding:dtDecoding:length(stim)*dtDecoding;
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

dt = dtDecoding*.001; %converts dt into seconds

numCells = size(kMat,2);
slen = size(histTermMat,1);

neglogli = 0;
grad = zeros(slen, 1);
hessian = zeros(slen, slen);

for i = 1:numCells
    
    currSpikeVec = spikeVecMat(:,i);
    kVec = kMat(:,i);
    kVec(length(kVec)+1:slen+ length(kVec)) = zeros(slen,1);
    cellKMat = triu(toeplitz(kVec));
    cellKMat = cellKMat(1:slen, size(kMat,1):slen+size(kMat,1)-1);
    
    stimCompInput = cellKMat*stim;
    fullFxnInput = stimCompInput + histTermMat(:,i);
    
    [condInt dcondInt ddcondInt logCondInt dlogCondInt ddlogCondInt] = expfunAndLog(fullFxnInput);
    
    cellNeglogli = -currSpikeVec'*logCondInt + sum(condInt)*dt;
    cellGrad = cellKMat'*(dcondInt*dt - currSpikeVec.*dlogCondInt);
    cellHessian = cellKMat*diag(ddcondInt*dt - currSpikeVec.*ddlogCondInt)*cellKMat';
%     
    neglogli = neglogli + cellNeglogli;
    grad = grad + cellGrad;
    hessian = hessian + cellHessian;
end

if nargin > 3
    cd('C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\decoding');
    load('alphaStimCovMat');
    stimCovMatInv = inv(stimCovMat(1:slen, 1:slen));
    priorNeglogli = .5*(stim'*stimCovMatInv*stim);% + (slen/2)*log(2*pi) + .5*log(stimCovMatDet);
    priorGrad = stim'*stimCovMatInv;
    priorHessian = stimCovMatInv;
    
    neglogli = neglogli + priorNeglogli;
    grad = grad + priorGrad';
    hessian = hessian + priorHessian;
end

logli = -neglogli;
grad = -grad;
hessian = -hessian;
