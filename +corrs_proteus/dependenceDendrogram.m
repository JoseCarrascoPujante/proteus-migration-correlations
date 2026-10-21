function Z = dependenceDendrogram(results)
%dependenceDendrogram  Average-linkage tree (UPGMA, Sokal & Michener 1958; MATLAB linkage) of the
%   eleven metrics on d = 1 - max(R*_n, 0), with R*_n the bias-corrected squared distance correlation of
%   the pair averaged over the three scenarios. Writes data/derived/dependenceDendrogram.csv with one row per
%   merge, as linkage returns it: the two clusters joined (metrics 1 to 11 in the order of metricNames, and
%   the cluster formed at merge k numbered 11 + k) and the merge height.
    arguments
        results (1, 1) struct
    end
    derivedDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'derived');
    conds = fieldnames(results);
    nM = size(results.(conds{1}), 2);
    pairs = nchoosek(1:nM, 2);
    meanR = nan(nM);
    for p = 1:size(pairs, 1)
        v = 0;
        for c = 1:numel(conds)
            v = v + corrs_proteus.distance_correlation(results.(conds{c})(:, pairs(p, 1)), ...
                results.(conds{c})(:, pairs(p, 2)));
        end
        v = v / numel(conds);
        meanR(pairs(p, 1), pairs(p, 2)) = v;
        meanR(pairs(p, 2), pairs(p, 1)) = v;
    end
    d = 1 - max(meanR, 0);
    d(1:nM + 1:end) = 0;
    Z = linkage(squareform(d, 'tovector'), 'average');
    writetable(array2table(Z, 'VariableNames', {'Cluster1', 'Cluster2', 'Height'}), ...
        fullfile(derivedDir, 'dependenceDendrogram.csv'));
end
