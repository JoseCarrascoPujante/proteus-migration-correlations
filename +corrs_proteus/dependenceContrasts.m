function T = dependenceContrasts()
%dependenceContrasts  Between-scenario difference in the dependence of every non-definitional metric pair:
%   R*_n in one scenario minus R*_n in another, tested by permuting the scenario labels of the pooled cells
%   after z-scoring each metric within its scenario (Szekely & Rizzo 2013, doi:10.1016/j.jmva.2013.02.012;
%   Phipson & Smyth 2010, doi:10.2202/1544-6115.1585; Benjamini & Hochberg 1995,
%   doi:10.1111/j.2517-6161.1995.tb02031.x).
    nPerm = 9999;
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    D = load(fullfile(projectRoot, 'data', 'derived', 'results.mat'));
    conds = ["NoStimuli", "Galvanotaxis", "Chemotaxis"];
    sfx = ["NoStim", "Galv", "Chem"];
    names = string(corrs_proteus.metricNames());
    pairs = nchoosek(1:numel(names), 2);
    pairs = pairs(~ismember(pairs, [6 7; 6 8], 'rows'), :);
    contrasts = [1 2; 1 3; 2 3];
    % R*_n is invariant to the location and scale of each metric, so z-scoring within the scenario
    % leaves it unchanged and keeps scenario shifts in the metrics out of the permuted samples.
    Z = arrayfun(@(c) zscore(D.results.(c)), conds, 'UniformOutput', false);
    n = size(Z{1}, 1);
    rs = RandStream('twister', 'Seed', 4);
    relabel = zeros(nPerm, 2 * n);
    for k = 1:nPerm, relabel(k, :) = randperm(rs, 2 * n); end

    nRow = size(pairs, 1) * size(contrasts, 1);
    [Pair, Contrast] = deal(strings(nRow, 1));
    [rFirst, rSecond, p] = deal(zeros(nRow, 1));
    row = 0;
    for ip = 1:size(pairs, 1)
        for ic = 1:size(contrasts, 1)
            row = row + 1;
            a = contrasts(ic, 1); b = contrasts(ic, 2);
            x = [Z{a}(:, pairs(ip, 1)); Z{b}(:, pairs(ip, 1))];
            y = [Z{a}(:, pairs(ip, 2)); Z{b}(:, pairs(ip, 2))];
            rFirst(row) = corrs_proteus.distance_correlation(x(1:n), y(1:n));
            rSecond(row) = corrs_proteus.distance_correlation(x(n+1:end), y(n+1:end));
            observed = abs(rFirst(row) - rSecond(row));
            exceed = false(nPerm, 1);
            parfor k = 1:nPerm
                g = relabel(k, :);
                d = corrs_proteus.distance_correlation(x(g(1:n)), y(g(1:n))) ...
                    - corrs_proteus.distance_correlation(x(g(n+1:end)), y(g(n+1:end)));
                exceed(k) = abs(d) >= observed;
            end
            p(row) = (1 + sum(exceed)) / (nPerm + 1);
            Pair(row) = names(pairs(ip, 1)) + " vs " + names(pairs(ip, 2));
            Contrast(row) = sfx(a) + "-vs-" + sfx(b);
        end
    end
    T = table(Pair, Contrast, rFirst, rSecond, rFirst - rSecond, p, bhfdr(p), ...
        'VariableNames', {'Pair', 'Contrast', 'bcDcor_first', 'bcDcor_second', 'Difference', ...
        'p_perm', 'q_FDR'});
    T = sortrows(T, {'q_FDR', 'p_perm'});
    for vn = T.Properties.VariableNames
        if isnumeric(T.(vn{1})), T.(vn{1}) = corrs_proteus.roundSigMin3(T.(vn{1})); end
    end
    writetable(T, fullfile(projectRoot, 'data', 'derived', 'dependenceContrasts.csv'));
end

function q = bhfdr(p)
    p = p(:); n = numel(p); [ps, ix] = sort(p);
    adj = ps .* n ./ (1:n)'; adj = flipud(cummin(flipud(adj)));
    q = ones(n, 1); q(ix) = min(adj, 1);
end
