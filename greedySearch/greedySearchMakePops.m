function [popNew, newPopList] = greedySearchMakePops(popOld, neuronsToSim, currSearchIter)


newPopList = [];
if currSearchIter == 1
    newPopList = neuronsToSim;
    
else
    % for every pop in popOld, make new pops by taking it and adding all
    % neurons in neuronsToSim
    for j = 1:length(popOld)
        oldPopList(:,j) = popOld(j).dat;
    end
    cnt = 1;
    for j = 1:size(oldPopList,2)
        for i = 1:length(neuronsToSim)
            newPopList(:,cnt) = sort([oldPopList(:,j); neuronsToSim(i)]);
            cnt = cnt + 1;
        end
    end
end
% remove repeats from newPopList
newPopList = unique(newPopList.','rows').';
% create the popNew struct
for i = 1:size(newPopList,2)
    popNew(i).dat = newPopList(:,i);
end