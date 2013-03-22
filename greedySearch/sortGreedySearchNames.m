function [sortedFiles] = sortGreedySearchNames(fileLists)


nameBase = 'greedySearch';
typeBase = '.mat';

num = [];
for i = 1:length(fileLists)
    temp = strrep(fileLists(i).name, nameBase, []);
    num(i) = str2num(strrep(temp, typeBase, []));
end

[vals, inds] = sort(num);
sortedFiles = fileLists(inds);