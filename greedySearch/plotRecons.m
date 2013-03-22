%% plot the reconstructions corresponding to the optimal population for a given stimulus

stimInd = 35;

reconIter = zeros(slen, numSearchIters);
for i  = 1:numSearchIters
    popInd = bestPopIdAll(i).dat;
    reconInd = (popInd - 1) * numSignals + stimInd;
    reconIter(:,i) = reconsAll(i).out(reconInd).optStim;
end

colorCode = colormap(redblue(20));

figure;
plot(1:slen, outSignal(:,stimInd), 'k', 'LineWidth', 2);
hold on;
plotInds = [1 5 10];
colorInds = [1 15 20];
for i = 1:length(plotInds);
plot(reconIter(:,plotInds(i)), 'LineWidth', 2, 'Color', colorCode(colorInds(i),:));
end

axis tight;