
b = sort(bestPop, 1, 'ascend');
tempPopList = b(:);
[stimFiltDiffOut histFiltDiffOut biasDiffOut] = glmDiff(cellModels(tempPopList), 2);

stimDiffMatAll = squareform(stimFiltDiffOut);

stimInds = reshape(1:100, 10, 10);
stimDiffMat = zeros(10, 10);
pValMat = zeros(10,10);
expectSelfVals = 45;
for i = 1:10
    aInds = stimInds(:,i);
    self_mat = stimDiffMatAll(aInds, bInds);
    sparsity_pattern = logical(triu(self_mat,1));
    selfVals = self_mat(sparsity_pattern);
    if numel(selfVals) < expectSelfVals
        selfVals(end+1:expectSelfVals) = 0;
    end
    for j = 1:10
        bInds = stimInds(:,j);
        totVals = stimDiffMatAll(aInds, bInds);
        totVals = totVals(:);
        if i == j
            totVals = selfVals;
        end
        [p, h] = ranksum(selfVals, totVals);
        pValMat(i,j) = p;
        val = mean(totVals);
        stimDiffMat(i,j) = val;
    end
end

figure;
imagesc(stimDiffMat); colorbar;
[r, v] = find(pValMat < .05);
hold on;
plot(v, r, '*')

