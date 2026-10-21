function [metricsTable, componentsTable, prTable] = principalComponents(results)
%principalComponents  Principal component analysis (MATLAB pca) of the eleven metrics, z-scored over
%   the 150 cells pooled over the scenarios, on their correlation matrix. Writes to data/derived/
%     pcaComponents.csv       variance explained by every component (%)
%     pcaMetrics.csv          loadings of every metric on the first three components, and its contribution to
%                             them: the squared loadings, which sum to one within a component, averaged with
%                             the variance each component explains as weight (%)
%     participationRatio.csv  PR = (sum of the eigenvalues)^2 / (sum of their squares) of the correlation
%                             matrix (Recanatesi et al. 2022, doi:10.1016/j.patter.2022.100555), pooled and
%                             within each scenario
    arguments
        results (1, 1) struct
    end
    derivedDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'derived');
    names = string(corrs_proteus.metricNames())';
    conds = fieldnames(results);
    nPC = 3;                                        % components summarized by the contribution

    pooled = cell2mat(cellfun(@(c) results.(c), conds, 'UniformOutput', false));
    [coeff, ~, ~, ~, explained] = pca(zscore(pooled));
    weights = explained(1:nPC) / sum(explained(1:nPC));
    contribution = (coeff(:, 1:nPC).^2) * weights * 100;
    metricsTable = table(names, coeff(:, 1), coeff(:, 2), coeff(:, 3), contribution, ...
        'VariableNames', {'Metric', 'PC1', 'PC2', 'PC3', 'ContributionPct'});
    componentsTable = table((1:numel(explained))', explained, 'VariableNames', {'Component', 'ExplainedPct'});

    Scope = ["Pooled"; string(conds)];
    ParticipationRatio = zeros(numel(Scope), 1);
    ev = eig(corr(zscore(pooled)));
    ParticipationRatio(1) = sum(ev)^2 / sum(ev.^2);
    for c = 1:numel(conds)
        ev = eig(corr(results.(conds{c})));
        ParticipationRatio(1 + c) = sum(ev)^2 / sum(ev.^2);
    end
    prTable = table(Scope, ParticipationRatio);

    writetable(componentsTable, fullfile(derivedDir, 'pcaComponents.csv'));
    writetable(metricsTable, fullfile(derivedDir, 'pcaMetrics.csv'));
    writetable(prTable, fullfile(derivedDir, 'participationRatio.csv'));
end
