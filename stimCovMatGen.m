%% generate stimCovMats for alpha conv noise
crossCorrLen = 100;
tauAlpha = 3;

crossCorrAll = zeros(3, crossCorrLen);
tempStim = getNoisyStim(5000000, 1, tauAlpha);
crossCorr1 = xcov(tempStim, crossCorrLen, 'coeff'); %


tauAlpha = 5;
tempStim = getNoisyStim(5000000, 1, tauAlpha);
crossCorr2 = xcov(tempStim, crossCorrLen, 'coeff'); %

tauAlpha = 10;
tempStim = getNoisyStim(5000000, 1, tauAlpha);
crossCorr3 = xcov(tempStim, crossCorrLen, 'coeff'); %

% plot(1:201, [crossCorr crossCorr2 crossCorr3]')

crossCorrAll(1,1:30) = crossCorr1(101:130);
crossCorrAll(2,1:41) = crossCorr2(101:141);
crossCorrAll(3,1:78) = crossCorr3(101:178);

%% playing with adding a cauchy-stylized stim corr fxn

temp = cauchypdf(0:1:199, 0, 3.5);
temp = temp./max(temp);
clear crossCorrAll
newCrossCorrAll = zeros(3, 200);
newCrossCorrAll(:,1:100) = crossCorrAll;
crossCorrAll(2,1:200) = temp;

%% adding a natural stimulus
timePts = 0:.001:5;

%generate a fast stimulus corrVec
tempStim = getNoisyStim(50000, 1, tauAlpha);
corrVecFast = xcov(tempStim, crossCorrLen, 'coeff'); %
corrVecFast = corrVecFast(101:147);
corrVecFast(length(corrVecFast)+1:length(timePts)) = 0;



oscFreq = 8;
corrVecOsc = cos(timePts*2*pi*oscFreq).*exp(timePts/-200);

slen = length(timePts);
stimCovMat = toeplitz(corrVecOsc);

% stimCovMat = eye(slen+slenBurn);
numSignals= 20;
dtStim = 1;

cholMat = chol(stimCovMat);
outSignal = zeros(slen, numSignals);
for j = 1:numSignals
    outSignal(:,j) = cholMat*randn(slen,1);
end

outSignalOsc = outSignal;

stimCovMat = toeplitz(corrVecFast);
cholMat = chol(stimCovMat);
outSignal = zeros(slen, numSignals);
for j = 1:numSignals
    outSignal(:,j) = cholMat*randn(slen,1);
end

outSignalFast = outSignal;

validInds = 100:slen;
stim = outSignalOsc(validInds,:) + .6*outSignalFast(validInds,:);

stim = zscore(stim);

crossCorrLenMixed = 2000;

for i = 1:numSignals
    crossCorrMixedAll(:,i) = xcov(stim(:,i), crossCorrLenMixed, 'coeff'); %
end
crossCorrMixed = mean(crossCorrMixedAll,2);
crossCorrMixedVec = crossCorrMixed((crossCorrLenMixed + 1):end);

stimCovMat = toeplitz(crossCorrMixed((crossCorrLenMixed + 1):end));
cholMat = chol(stimCovMat);

outSignal = zeros(crossCorrLenMixed+1, numSignals);
for j = 1:numSignals
    outSignal(:,j) = cholMat*randn(crossCorrLenMixed+1,1);
end

outSignalMixed = outSignal;

validInds = 200:crossCorrLenMixed+1;
outSignalMixed = zscore(outSignalMixed(validInds,:));

[V,D]=eig(A);

d=diag(D);
d(d<=0)=eps;

Anew= V*diag(d)*V';

tapers=3; avg=1; fs = 1000;
pars = struct ('Fs', fs, ...
    'tapers', [tapers 2*tapers-1], ...
    'pad', 0, ...
    'trialave', 1);
[S, f] = mtspectrumc(stim, pars);

figure; plot(f,S); ylabel('power'); xlabel('Frequency (Hz)');
loglog(f,S)
plot(f,S)
%

% 
% outSignal = outSignal(slenBurn+1:end,:);
% stimCovMat= stimCovMat(slenBurn+1:end, slenBurn+1:end);

%% try making a stimCovMat directly from a signal with properties i like
 + cos([.001:.001:500]*2*pi*oscFreq)';
tempStim =  getNoisyStim(500000, 1, 10);
% plot(tempStim)
crossCorrNat = xcov(tempStim, 1000, 'coeff'); %
% crossCorrMixed = mean(crossCorrMixedAll,2);
crossCorrMixedVec = crossCorrNat(((length(crossCorrNat)-1)/2)+1:end);

crossCorrLenMixed = length(crossCorrMixedVec);

stimCovMat = toeplitz(crossCorrMixedVec);
cholMat = chol(stimCovMat);

numSignals = 20;
outSignal = zeros(crossCorrLenMixed, numSignals);
for j = 1:numSignals
    outSignal(:,j) = cholMat*randn(crossCorrLenMixed,1);
end

outSignalMixed = outSignal;

tapers=3; avg=1; fs = 1000;
pars = struct ('Fs', fs, ...
    'tapers', [tapers 2*tapers-1], ...
    'pad', 0, ...
    'trialave', 1);
[S, f] = mtspectrumc(outSignal, pars);

figure; plot(f,S); ylabel('power'); xlabel('Frequency (Hz)');
loglog(f,S)

%% make stimCovMat from 1/f dist'n using powernoise fxn

x = powernoise(1, 1000, 'normalize');

%% generate OU stim

gVals = [.05 .1 .2];
tau = [3 10];
gVals = 1./tau;

slen = 5000000;
dt = 1;
t = dt:dt:slen;

D = 10;

crossCorrLen = 1000;

for i = 1:length(gVals)
    g = gVals(i);
    
    X_free = cumsum(sqrt(2*D*dt)*randn(size(t)));
    X_filt = filter([0 g*dt], [1 -1+g*dt] , X_free);
    X_OU = zscore(X_free-X_filt);
    tempStim = X_OU;
    corrVecOU = xcov(tempStim, crossCorrLen, 'coeff'); %
    crossCorrAll(i,:) = corrVecOU(((length(corrVecOU)-1)/2)+1:end);
end

% or just generate the auto corr fxn exactly:
% exp(1/(-tau)); % tau = corr time
tau = [20 40];
crossCorrTPts = 0:1:1000;
for i = 1:length(tau)
    crossCorrAll(i,:) = exp(-crossCorrTPts/(tau(i)));
end

[outSignal, stimCovMat, stimParms] = generateStims(crossCorrAll, 512, 50);
%% generate oscillatory stim

timePts = 0:.001:5;

oscFreqs = [8 20 60];
for i = 1:length(oscFreqs)
    oscFreq = oscFreqs(i);
    crossCorrAll(i,:) = cos(timePts*2*pi*oscFreq).*exp(timePts/-1000);
end

%% regenerate a naturalistic stim
tau = 10;
timePts = 0:1:5000;
corrVecFast  = exp(-timePts/tau);

oscFreq = .008;
corrVecOsc = cos(timePts*2*pi*oscFreq).*exp(timePts/-200);

slen = length(timePts);
stimCovMat = toeplitz(corrVecOsc);

% stimCovMat = eye(slen+slenBurn);
numSignals= 20;
dtStim = 1;

cholMat = chol(stimCovMat);
outSignal = zeros(slen, numSignals);
for j = 1:numSignals
    outSignal(:,j) = cholMat*randn(slen,1);
end

outSignalOsc = outSignal;

stimCovMat = toeplitz(corrVecFast);
cholMat = chol(stimCovMat);
outSignal = zeros(slen, numSignals);
for j = 1:numSignals
    outSignal(:,j) = cholMat*randn(slen,1);
end

outSignalFast = outSignal;

validInds = 100:slen;
stim = outSignalOsc(validInds,:) + .6*outSignalFast(validInds,:);

stim = zscore(stim);

crossCorrLenMixed = 2000;

parfor i = 1:numSignals
    crossCorrMixedAll(:,i) = xcov(stim(:,i), crossCorrLenMixed, 'coeff'); %
end
crossCorrMixed = mean(crossCorrMixedAll,2);
crossCorrMixedVec = crossCorrMixed((crossCorrLenMixed + 1):end);

stimCovMat = toeplitz(crossCorrMixedVec);
cholMat = chol(stimCovMat);

outSignal = zeros(crossCorrLenMixed+1, numSignals);
for j = 1:numSignals
    outSignal(:,j) = cholMat*randn(crossCorrLenMixed+1,1);
end

outSignalMixed = outSignal;

validInds = 100:crossCorrLenMixed+1;
outSignalMixed = zscore(outSignalMixed(validInds,:));

[V,D]=eig(A);

d=diag(D);
d(d<=0)=eps;

Anew= V*diag(d)*V';

tapers=3; avg=1; fs = 1000;
pars = struct ('Fs', fs, ...
    'tapers', [tapers 2*tapers-1], ...
    'pad', 0, ...
    'trialave', 1);
[S, f] = mtspectrumc(stim, pars);

figure; plot(f,S); ylabel('power'); xlabel('Frequency (Hz)');
loglog(f,S)
plot(f,S)
%

%% another natural stim
osiLen = 10000;
rate = 50;
tau = 10;
oscFreq = 8;
epspAmp = 25;
oscAmp = 50;

clear stim
numSignals = 50;

for i = 1:numSignals
[outputCurrent] = osiStimGen2(osiLen, rate, tau, oscFreq, epspAmp, oscAmp);
outputCurrent = downsample(outputCurrent, 10);
stim(:,i) = outputCurrent;
end
validInds = 100:osiLen;
stim = stim(validInds,:);
% stim = stim + 3*filter(Hd2, randn(length(validInds), numSignals));
stim = zscore(stim);



clear crossCorrMixedAll
crossCorrLenMixed = 1000;

parfor i = 1:numSignals
    crossCorrMixedAll(:,i) = xcov(stim(:,i), crossCorrLenMixed, 'coeff'); %
end
crossCorrMixed = mean(crossCorrMixedAll,2);
crossCorrMixedVec = crossCorrMixed((crossCorrLenMixed + 1):end);

stimCovMat = toeplitz(crossCorrMixedVec);
cholMat = chol(stimCovMat);

outSignal = zeros(crossCorrLenMixed+1, numSignals);
for j = 1:numSignals
    outSignal(:,j) = cholMat*randn(crossCorrLenMixed+1,1);
end

outSignalMixed = outSignal;

validInds = 100:crossCorrLenMixed+1;
outSignalMixed = zscore(outSignalMixed(validInds,:));
