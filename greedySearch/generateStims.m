function [outSignal, stimCovMat, stimParms] = generateStims(crossCorrFxn, stimLen, numSignalsPerStatIn, varargin)

if nargin > 3
    if varargin{1} == 1
        %     vargin{1} = randSeed;
        s = RandStream('mcg16807', 'Seed',0);
        RandStream.setDefaultStream(s);
    end
end

% load 'C:\Users\Shreejoy\Desktop\shreejoy\neural_diversity\ST_machinelearn\reconKP\populationsims\oscStimCovMat';
% crossCorrAll = crossCorrAll(2,:);
crossCorrAll = crossCorrFxn;
slen = stimLen;

slenBurn = 500;
numSignalsPerStat = numSignalsPerStatIn;
numStimStats = size(crossCorrAll,1);

numSignals= numSignalsPerStat*numStimStats;
dtStim = 1;


stimCovMat = zeros(slen+slenBurn, slen+slenBurn, numStimStats);
outSignal = zeros(slen+slenBurn, numSignals);
for k = 1:numStimStats
    corrVec = nonzeros(crossCorrAll(k,:))';
    
    if length(corrVec) < slen+slenBurn
        corrVec(end+1:slen+slenBurn) = 0;
    else
        corrVec = corrVec(1:slen+slenBurn);
    end
    stimCovMat(:,:,k) = toeplitz(corrVec);
    
    
    cholMat = chol(stimCovMat(:,:,k));
    
    % not supporting dt's that aren't 1 ms yet
    for j = 1:numSignalsPerStat
        outSignal(:,j+(k-1)*numSignalsPerStat) = cholMat*randn(slen+slenBurn,1);
    end
    
end
% outSignal = outSignal(slenBurn+1:end,:);
% stimCovMat= stimCovMat(slenBurn+1:end, slenBurn+1:end,:);

validInds = slenBurn/2+1 :slen+slenBurn -slenBurn/2;
outSignal = outSignal(validInds,:);
stimCovMat= stimCovMat(validInds, validInds,:);

stimParms.outSignal = outSignal;
stimParms.numSignals = numSignals;
stimParms.numStimStats = numStimStats;

for i = 1:numStimStats
    stimEnt(i) = .5*sum(myLog2(eig(stimCovMat(:,:,i)))) + .5*slen*myLog2(2*pi*exp(1));
end
% % 
% tapers=3; avg=1; fs = 1000;
% pars = struct ('Fs', fs, ...
%     'tapers', [tapers 2*tapers-1], ...
%     'pad', 0, ...
%     'trialave', 1);
% for k = 1:numStimStats
%     [S(:,k), f] = mtspectrumc(outSignal(:,1+(k-1)*numSignalsPerStat:k*numSignalsPerStat), pars);
% end
% figure; plot(f,S); ylabel('power'); xlabel('Frequency (Hz)');
% loglog(f,S); 
% plot(f,S)
