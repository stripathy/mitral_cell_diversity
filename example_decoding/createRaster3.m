function [] = createRaster3(cellStats, preStimTime, postStimTime, colorVals)
%scavenged from LstStitchEnsRaster02 from EST at DalyLab summer 2008
% plots raster in current figure window
% colorVals is a numTrials x 
% same as 

lineWidth = 2;
APTraces = cellStats;
if isstruct(APTraces)
    %--- PSTH with TIME stamps:
    ctr=1;
    cumLen = length(cellStats);
    %---
    % FigPos=[60 10 500 400]; % [left top width height]
    % figure('Position',FigPos)
    % set(gcf,'Color','w');
    % h(1)=axes('Position',[.11 .08 .85 .85]); % RECT = [left, bottom, width, height]
    hold on
    for ctr=1:cumLen
        if isempty(APTraces(:,ctr).list) % This is to ACCOUNT for EMPTY MATRIX cases!!!
            if ctr == cumLen
                CurrAPTs(:,1)=preStimTime;
                CurrAPn = 1;
                CurrCtrA(1:CurrAPn,1)=ctr;
                CtrA=1;
                for CtrA=1:CurrAPn
                    line([CurrAPTs(CtrA,1) CurrAPTs(CtrA,1)],[CurrCtrA(CtrA,1)-.35 CurrCtrA(CtrA,1)+.35],'Color','b','LineWidth',.5)
                end
                clear CurrCtrA CurrAPn CurrAPTs
            end
            ctr = ctr + 1;
        else
            CurrAPTs(:,1)=APTraces(ctr).list(find(APTraces(ctr).list >= preStimTime & APTraces(ctr).list < postStimTime));
            %         CurrAPTs(:,1)=APTraces(ctr).spkVector(find(APTraces(ctr).spkVector > preStimTime & APTraces(ctr).spkVector < postStimTime),ctr);
            CurrAPn=length(CurrAPTs);
            CurrCtrA(1:CurrAPn,1)=ctr;
            CtrA=1;
            for CtrA=1:CurrAPn
                line([CurrAPTs(CtrA,1) CurrAPTs(CtrA,1)],[CurrCtrA(CtrA,1)-.35 CurrCtrA(CtrA,1)+.35],'Color','b','LineWidth',.5)
            end
            clear CurrCtrA CurrAPn CurrAPTs
        end
    end
else
    APTraces = APTraces';
        %--- PSTH with TIME stamps:
    ctr=1;
    [cumLen blah] = size(cellStats);
    %---
    % FigPos=[60 10 500 400]; % [left top width height]
    % figure('Position',FigPos)
    % set(gcf,'Color','w');
    % h(1)=axes('Position',[.11 .08 .85 .85]); % RECT = [left, bottom, width, height]
    hold on
    for ctr=1:cumLen
        if nnz(APTraces(:,ctr)) == 0 % This is to ACCOUNT for EMPTY MATRIX cases!!!
            if ctr == cumLen
                CurrAPTs(:,1)=preStimTime;
                CurrAPn = 1;
                CurrCtrA(1:CurrAPn,1)=ctr;
                CtrA=1;
                for CtrA=1:CurrAPn
                    line([CurrAPTs(CtrA,1) CurrAPTs(CtrA,1)],[CurrCtrA(CtrA,1)-.35 CurrCtrA(CtrA,1)+.35],'Color',colorVals(ctr,:),'LineWidth',lineWidth)
                end
                clear CurrCtrA CurrAPn CurrAPTs
            end
            ctr = ctr + 1;
        else
            temp = nonzeros(APTraces(:,ctr));
            inds = temp >= preStimTime & temp < postStimTime;
            CurrAPTs(:,1)=APTraces(inds,ctr);
            %         CurrAPTs(:,1)=APTraces(ctr).spkVector(find(APTraces(ctr).spkVector > preStimTime & APTraces(ctr).spkVector < postStimTime),ctr);
            CurrAPn=length(CurrAPTs);
            CurrCtrA(1:CurrAPn,1)=ctr;
            CtrA=1;
            for CtrA=1:CurrAPn
                line([CurrAPTs(CtrA,1) CurrAPTs(CtrA,1)],[CurrCtrA(CtrA,1)-.35 CurrCtrA(CtrA,1)+.35],'Color',colorVals(ctr,:),'LineWidth',lineWidth)
            end
            clear CurrCtrA CurrAPn CurrAPTs
        end
    end
end
set(gca,'YDir','reverse');
axis([-preStimTime postStimTime 0 cumLen+1 ]);