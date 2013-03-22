%decoding test on GWN stim from LNP neurons



%% 1.  Set parameters and display for GLM % =============================
RefreshRate = 1000; 
DTsim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 50;    % Number of time bins in stimulus filter
ttk = [-nkt+1:0]';  % time relative to spike of stim filter taps
ggsim = makeSimStruct_GLM(nkt,DTsim); % Creat
ggsim.RefreshRate = RefreshRate;

swid = 1;
slen = 100; % Stimulus length (frames);  More samples gives better fit
% Stim = round(rand(slen,swid))*2-1;  %  Run model on long, binary stimulus
Stim = randn(slen,swid);
[tsp,vmem,ispk] = simGLM(ggsim,Stim);  % run model
nsp = length(tsp);

stimInterval = [1 slen];
dtStim = 1;
stimCovMat = eye(slen);

gg = ggsim;

cellModels = gg;
spikeTimes = tsp;
[optStim, exitflag, fval] = bayesStimDecoder1(cellModels, spikeTimes, stimInterval, dtStim, stimCovMat);

