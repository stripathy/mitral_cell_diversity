function [logli, grad, hessian] = bayesDecoder_logli(prs, inputParms)
% computes the most likely stimulus given spike trains from a number of
% cells, an encoding model for the cells, and a prior of the stim distn

% inputs in some fashion:
% 
% vectorized spikes of ith cell during stim length
% distn over stimuli
% stim length, dt of stim
% nonlinearity fxn f
% encoding model for each cell using to decode
spikeMat = inputParms.spikeMat; %rows are spikes in bin i, cols are diff cells
numCells = size(spikeMat,2);
slen = size(spikeMat,1);
currEstStim = prs;
dtStim = inputParms.dtStim*.001;

kFiltMat = inputParms.kFiltMat; % matrix where filters of ith cell are along column
staWin = size(kFiltMat, 1);

%actually I should just send a matrix of length stim with these values in
%it plus a vector for the offsets
hFiltFact = inputParms.hFiltFact;% matrix where filters of ith cell are along column

% priorMeans = zeros(1, slen); %vector of means (should be zeros for zero-mean stim
priorCovMatDet = inputParms.priorCovMatDet;
priorInvCovMat = inputParms.priorInvCovMat;

logli = 0;
grad = zeros(slen, 1);
hessian = zeros(slen,slen);
% I can actually precompute the stuff not dependent on the current stim
% compute logli and sum for each cell (maybe i can use vectorized ops x
% cell)    
for i = 1:numCells
    currKFilt = kFiltMat(:,i);
    fullKFiltVec = zeros(1, slen+staWin);
%     fullKFiltVec = zeros(1, slen);
    fullKFiltVec(1:staWin) = kFiltMat(:,i);
    kMat = toeplitz(fullKFiltVec);
    
    kMat = kMat(1:slen, staWin:slen+staWin-1);
    kMat = triu(kMat, -(staWin-1));
    
    % evaluate the input to the function f for each cell
    fxnValStim = kMat * currEstStim;
    fxnInput = fxnValStim+hFiltFact(:,i);
    
%     plot(1:250, fxnValStim, 1:250, hFiltFact)
    
    [fxnVal fxnValDer fxnValDerDer logfxnVal logfxnValDer logfxnValDerDer] =...
        expfunAndLog(fxnInput); %hard coding in expfun
    
    cellSpikeVec = spikeMat(:,i); % make sure to get the time offset just right!
    cellLogli = -cellSpikeVec'*logfxnVal + sum(fxnVal*dtStim);
    cellGrad = kMat*(fxnValDer*dtStim - cellSpikeVec.*logfxnValDer);
    cellHessian = kMat'*diag(fxnValDerDer*dtStim-cellSpikeVec.*logfxnValDerDer)*kMat;
    
    if sum((fxnValDerDer*dtStim-cellSpikeVec.*logfxnValDerDer)<0)>1
        a = 1;
    end
    
    logli = logli + cellLogli;
    grad = grad + cellGrad;
    hessian = hessian + cellHessian;
end

% compute logli for prior stim distn
priorLogli = (.5*slen*log(2*pi) + .5*log(priorCovMatDet)+ .5*currEstStim'*priorInvCovMat*currEstStim); %already negated
priorGrad = currEstStim'*priorInvCovMat;
priorHessian = priorInvCovMat';
%     
logli = logli + priorLogli;
grad = grad + priorGrad';
hessian = hessian + priorHessian;

    