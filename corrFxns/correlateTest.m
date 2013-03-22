
% spikes should be a matix with a row for each cell in the network,
% spike times are listed in the columns.  multiple trials or realizations
% can be listed in 3rd axis

function [rhos,acorrs,xcorrs,Instrate] = correlateTest(spikes,tmax,trialsPerCell)

span=25;  % length (in ms) on either side of the correlograms
pad=0;     % size of box fxn added to each side of spike
T=4;       % span (ms) over which area under correlograms is integrated

numCells=size(spikes,1);  % # neurons
reps = trialsPerCell;
% reps=size(spikes,3); % # trials
rhos=zeros(numCells,numCells);    % matrix of pair-wise spike correlations
acorrs=zeros(numCells,span*2+1);  % matrix of all autocorrelograms
xcorrs=zeros(numCells,span*2+1,numCells); % matrix of all cross-correlograms
Instrate=zeros(numCells,tmax,max(reps)); % binary rep. of rate
numSpikes=zeros(numCells,max(reps)); %rates for each cell, each trial
firingRates = zeros(numCells,1);

%pad the spikes and turn to instantaneous rate
for aa=1:numCells
    for dd=1:reps(aa)
        tempNumSpikes=nnz(spikes(aa,:,dd));
        numSpikes(aa,dd)=tempNumSpikes;
        if tempNumSpikes>0
            k=int32(round(spikes(aa,1:tempNumSpikes,dd)));
            if k(1) == 1 || k(1) == 0
                k(1) = 1;
            end
            Instrate(aa,k,dd)=1;
            if pad>0
                for ll=1:pad
                    Instrate(aa,k-ll,dd)=1;
                    Instrate(aa,k+ll,dd)=1;
                end
            end
        end
    end
    %compute firingRate
    firingRates(aa) = mean(numSpikes(aa,:))/(tmax*.001);
end  %makes matrix of instantaneous rate (binary) to be used for correlation function


%compares spike times to generate all pairwise correlograms
for bb=1:numCells
    for ee=1:reps(bb) %num repeats cell 1
        for cc=1:numCells
            for ff=1:reps(cc) %num repeats cell 2
                for dd=1:numSpikes(bb,ee) %for all spikes in cell 1
                    currSpikeTime=round(spikes(bb,dd,ee)); %find every spike in bbth cell
                    if currSpikeTime>span && currSpikeTime<(tmax-span)
                        xcorrs(bb,:,cc)=xcorrs(bb,:,cc)+Instrate(cc,currSpikeTime-span:currSpikeTime+span,ff)/(reps(bb)*reps(cc)); %/numSpikes(cc,ff)...
                            %-nnz(Instrate(cc,currSpikeTime-span:currSpikeTime+span,ff))/((2*span+1)*numSpikes(cc,ff)); %check out this subtracting part
                        %subtract out count1*binsize*rate2
                    end
                end
            end
        end
    end
end



%     xcorrs=xcorrs./reps;
for aa=1:numCells
    for bb = 1:numCells
        xcorrs(aa,:,bb)= xcorrs(aa,:,bb)/.001 - sqrt(firingRates(aa)*firingRates(bb));
        xcorrs(aa,:,bb)=xcorrs(aa,:,bb)/(sqrt((firingRates(aa)*firingRates(bb)))*.001);
    end
    acorrs(aa,:)=xcorrs(aa,:,aa);
end
    
    %integrates curves to calculate correlation coeff
    
    areaX=zeros(numCells,numCells);
    areaA=zeros(numCells,1);

    for ee=1:numCells
        for ff=1:numCells
            for dd=1:(2*T+1)
                if ee==ff
                    areaA(ee,1)=areaA(ee,1)+acorrs(ee,span-T+dd)*((T-abs(-T+dd))/T);
                end
                    areaX(ee,ff)=areaX(ee,ff)+xcorrs(ee,span-T+dd,ff)*((T-abs(-T+dd))/T);
            end
        end
    end
    
    for aa=1:numCells
        for bb=1:numCells
            rhos(aa,bb)=areaX(aa,bb)/(sqrt(areaA(aa)*areaA(bb)));
        end
    end
    
    
end