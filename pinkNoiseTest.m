x = spatialPattern([1000 1],-1);
x = zscore(x);

xResamp = fastinterp2(x, 1000);

xResamp = resample(x, 1000, 1);
xResamp = zscore(xResamp);

plot(xResamp)

[pxx, f] = pwelch(xResamp, [], [], [], 1000);
loglog(f, pxx)

[pxx, f] = pwelch(x, [], [], [], 1000);
loglog(f, pxx, [0 500], [0 10^-(log10(500))])

plot(f, pxx)

sigLen = 1000;
dt = 1/1000;
plot(dt:dt:sigLen, xResamp)

whiteNoise = randn(10000000, 1);
windowSize = 10;
brownNoise = filter(ones(1,windowSize)/windowSize,1,whiteNoise);
brownNoise = zscore(brownNoise);
crossCorrLen = 500;
crossCorr = xcov(brownNoise, crossCorrLen, 'coeff'); 


% brownNoise = cumsum(randn(10000,1));
% brownNoise = zscore(brownNoise);
[pxx, f] = pwelch(brownNoise, [], [], [], 1000);
loglog(f, pxx)

crossCorrTemp = crossCorr(501:500+19);

corrVec = nonzeros(crossCorrTemp)';
    corrVec(end+1:1024) = 0;
    stimCovMat = toeplitz(corrVec);
    
    
    cholMat = chol(stimCovMat);
    genSignal = cholMat*randn(1024,1);
    % not supporting dt's that aren't 1 ms yet
    for j = 1:100
        genSignal(:,j) = cholMat*randn(1024,1);
    end
    tapers=3; avg=1; fs = 1000;
pars = struct ('Fs', fs, ...
    'tapers', [tapers 2*tapers-1], ...
    'pad', 0, ...
    'trialave', 1);
[S, f] = mtspectrumc(genSignal, pars);
loglog(f,S)

for k = 1:numSignals
    [S(:,k), f] = mtspectrumc(outSignal(:,k), pars);
end
loglog(f,S)
    