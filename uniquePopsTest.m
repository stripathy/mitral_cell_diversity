cellMat = zeros(11^5, 5);

cnt = 1;
for i = 1:11
    for j = 1:11
        for k = 1:11
            for m = 1:11
                for n = 1:11
                    cellMat(cnt,:) = [i j k m n];
                    cnt = cnt + 1;
                end
            end
        end
    end
end

cellMat = sort(cellMat,2);

newCellMat = unique(cellMat, 'rows');

% how many of these have more than 3 repeats?
b = zeros(size(newCellMat));
cnt = 1; clear validPops
for i = 1:size(newCellMat,1)
    temp = unique(newCellMat(i,:));
    if length(temp) > 3
        validPops(cnt) = i;
        cnt = cnt + 1;
    end
    b(i,1:length(temp)) = temp;
end

nonRepCellMat = newCellMat(validPops,:);

