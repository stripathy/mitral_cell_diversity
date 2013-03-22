function [filtSpikes] = getSpikesInt(spikes, timeInterval)

tempTestSpikes = zeros(4000,size(spikes,2));
for i = 1:size(spikes,2)
    temp = spikes(spikes(:,i)>timeInterval(1) & spikes(:,i) < timeInterval(2),i);
    tempTestSpikes(1:length(temp),i) = temp;
end

maxSpikes = 4000 - min(sum(tempTestSpikes==0));

filtSpikes(1:maxSpikes,:) = tempTestSpikes(1:maxSpikes,:);

if isempty(filtSpikes)
    filtSpikes = 0;
end