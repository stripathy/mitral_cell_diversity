% test decoding for gaussian approx to the likelihood

neurons = [10 28 33 44 44];
pop.dat = neurons;

% compute the mle's and hessians for each neuron of a pop independently
reconStructInd = 5006;
j = reconStructInd;
pop.dat = reconStructIn(j).popMakeup;
initStim = outSignal(:,mod(j-1,numSignals)+1);
testSpikes = reconStructIn(j).testSpikes;  
stimCovMatIn = -1;
stimInterval = [0 slen];
parfor i = 1:length(reconStructIn(j).popMakeup)
    
%     tempPop.dat = pop.dat(i);
    currPop = cellModels(pop.dat(i));
    tempSpikes = testSpikes(:,i);
    stimCovMatIn = -1;
%     
    [optStim, exitflag, likeli, hessian] = bayesStimDecoder1(currPop, tempSpikes, stimInterval, dtStim, stimCovMatIn, initStim);
    reconStructOut(i).optStim = optStim;
    reconStructOut(i).likeli = likeli;
    reconStructOut(i).hessian = hessian;
%     reconStructOut(i).stimErrorBars = sqrt(diag(inv(hessian)));
end  


sMapApprox
currStimCovMat = stimCovMat(:,:,1);
hessTerm = zeros(slen, slen);
stimTerm = zeros(slen,1);

for i = 1:length(neurons)
    tempHess = reconStructOut(i).hessian;
    tempStim = reconStructOut(i).optStim;
    stimTerm = stimTerm + tempHess*tempStim;
    hessTerm = hessTerm + tempHess;
end
sMapApprox =  (eye(slen) + currStimCovMat * hessTerm)\ (currStimCovMat * stimTerm);

sMapApprox = (inv(currStimCovMat) + hessTerm)\(stimTerm);

figure; 
plot(1:slen, initStim, 'k', 1:slen, sMapApprox, 'r');

sMap = (eye(slen) + stimCovMat(:,:,1)*hessian)\(stimCovMat(:,:,1)*hessian*optStim);