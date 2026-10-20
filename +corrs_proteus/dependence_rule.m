function T = dependence_rule(results)
%dependence_rule  Dependence of every metric pair: distance covariance permutation test and bias-corrected
%   squared distance correlation with BCa 95% CI (Szekely et al. 2007, doi:10.1214/009053607000000505;
%   Szekely & Rizzo 2013, doi:10.1016/j.jmva.2013.02.012; DiCiccio & Efron 1996, doi:10.1214/ss/1032280214;
%   permutation p-value as in Phipson & Smyth 2010, doi:10.2202/1544-6115.1585),
%   RESET form test (Ramsey 1969, doi:10.1111/j.2517-6161.1969.tb00796.x; Ramsey & Schmidt 1976,
%   doi:10.1080/01621459.1976.10480355), U-shape test (Lind & Mehlum 2010, doi:10.1111/j.1468-0084.2009.00569.x),
%   Benjamini-Hochberg FDR (doi:10.1111/j.2517-6161.1995.tb02031.x).
    arguments
        results (1, 1) struct
    end
    conds = string(fieldnames(results))';
    assert(isequal(conds, ["NoStimuli", "Galvanotaxis", "Chemotaxis"]), ...
        'dependence_rule:conditions', 'Expected the scenarios NoStimuli, Galvanotaxis, Chemotaxis.');
    sfx = ["NoStim", "Galv", "Chem"];
    dirs = ["AxB", "BxA"];
    names = string(corrs_proteus.metricNames());
    assert(all(structfun(@(M) all(isfinite(M), 'all'), results)), 'dependence_rule:missing', ...
        'results hold NaN or Inf: every metric must be defined for every cell.');
    nCond = numel(conds);
    pairs = nchoosek(1:numel(names), 2);
    nPair = size(pairs, 1);
    alpha = 0.05;

    isDef = ismember(pairs, [6 7; 6 8], 'rows');

    nBoot = 2000;     % bootstrap resamples of the BCa interval of R*_n
    bootOpts = statset('UseParallel', false, 'Streams', RandStream('twister', 'Seed', 2));
    dcorStat = @(x, y) corrs_proteus.distance_correlation(x, y);
    nPerm = 9999;
    rsPerm = RandStream('twister', 'Seed', 1);
    [bcDc, pDc, bcLo, bcHi, pe, ppe, sp, psp] = deal(zeros(nPair, nCond));
    for c = 1:nCond
        M = results.(conds(c));
        for p = 1:nPair
            a = M(:, pairs(p, 1));
            b = M(:, pairs(p, 2));
            bcDc(p, c) = corrs_proteus.distance_correlation(a, b);
            pDc(p, c) = dcovPermutationTest(a, b, nPerm, rsPerm);
            ci = bootci(nBoot, {dcorStat, a, b}, 'Type', 'bca', 'Alpha', alpha, 'Options', bootOpts);
            bcLo(p, c) = ci(1);
            bcHi(p, c) = ci(2);
            [pe(p, c), ppe(p, c)] = corr(a, b, 'Type', 'Pearson');
            % Spearman rho as the Pearson correlation of the tied ranks, tested with its t statistic:
            % MATLAB's 'Spearman' p-value approximation returns 0 far in the tail.
            [sp(p, c), psp(p, c)] = corr(tiedrank(a), tiedrank(b));
        end
    end
    qDc = reshape(bhfdr(pDc(:)), nPair, nCond);
    qPe = reshape(bhfdr(ppe(:)), nPair, nCond);
    qSp = reshape(bhfdr(psp(:)), nPair, nCond);
    nCondDc = sum(qDc < alpha, 2);
    dependent = nCondDc == nCond;
    coupled = dependent & ~isDef;

    % Form and trend-reversal tests for every pair, in both regression directions, corrected as one
    % family of 330 tests so that every pair is judged on the same footing.
    [pReset, pU] = deal(NaN(nPair, nCond, 2));
    for k = 1:nPair
        for c = 1:nCond
            M = results.(conds(c));
            a = M(:, pairs(k, 1));
            b = M(:, pairs(k, 2));
            [pReset(k, c, 1), pU(k, c, 1)] = formTests(a, b);
            [pReset(k, c, 2), pU(k, c, 2)] = formTests(b, a);
        end
    end
    qReset = reshape(bhfdr(pReset(:)), size(pReset));
    qU = reshape(bhfdr(pU(:)), size(pU));

    form = repmat("linear", nPair, 1);
    form(any(qReset < alpha, [2 3])) = "nonlinear";
    nReversal = sum(qU < alpha, [2 3]);
    pairRank = zeros(nPair, 1);
    [~, order] = sort(median(bcDc, 2), 'descend');
    pairRank(order) = (1:nPair)';

    pairName = names(pairs(:, 1))' + " vs " + names(pairs(:, 2))';
    T = table(pairName, isDef, dependent, coupled, nCondDc, pairRank, form, nReversal, ...
        median(bcDc, 2), ...
        'VariableNames', {'Pair', 'Definitional', 'Dependent', 'Coupled', 'nCond_dCov', ...
        'Rank', 'Form', 'nReversal', 'bcDcor_median'});
    for c = 1:nCond
        s = sfx(c);
        T.(char("bcDcor_" + s)) = bcDc(:, c);
        T.(char("bcDcor_lo_" + s)) = bcLo(:, c);
        T.(char("bcDcor_hi_" + s)) = bcHi(:, c);
        T.(char("pDcov_" + s)) = pDc(:, c);
        T.(char("qDcov_" + s)) = qDc(:, c);
        T.(char("Pearson_" + s)) = pe(:, c);
        T.(char("qPearson_" + s)) = qPe(:, c);
        T.(char("Spearman_" + s)) = sp(:, c);
        T.(char("qSpearman_" + s)) = qSp(:, c);
        for d = 1:2
            T.(char("qRESET_" + dirs(d) + "_" + s)) = qReset(:, c, d);
        end
        for d = 1:2
            T.(char("qLindMehlum_" + dirs(d) + "_" + s)) = qU(:, c, d);
        end
    end
end

function [pReset, pU] = formTests(x, y)
    xs = (x - mean(x)) / std(x);
    ys = (y - mean(y)) / std(y);
    n = numel(ys);
    X0 = [ones(n, 1) xs];
    yHat = X0 * (X0 \ ys);
    sse0 = sum((ys - yHat) .^ 2);
    X1 = [X0 yHat .^ 2 yHat .^ 3 yHat .^ 4];
    nQ = 3;
    sse1 = sum((ys - X1 * (X1 \ ys)) .^ 2);
    df1 = n - size(X1, 2);
    pReset = fcdf(((sse0 - sse1) / nQ) / (sse1 / df1), nQ, df1, 'upper');

    X2 = [X0 xs .^ 2];
    b2 = X2 \ ys;
    df2 = n - 3;
    covB = (sum((ys - X2 * b2) .^ 2) / df2) * ((X2' * X2) \ eye(3));
    ends = [min(xs) max(xs)];
    slope = b2(2) + 2 * b2(3) * ends;
    se = sqrt(covB(2, 2) + 4 * ends * covB(2, 3) + 4 * ends .^ 2 * covB(3, 3));
    tEnds = slope ./ se;
    s = sign(b2(3));
    pU = max(tcdf(s * tEnds(1), df2), tcdf(-s * tEnds(2), df2));
end

function pValue = dcovPermutationTest(x, y, nPerm, rs)
    % Distance covariance permutation test: the observed V-statistic against nPerm permutations of y.
    n = numel(x);
    a = abs(x - x.');
    b = abs(y - y.');
    A = a - mean(a, 2) - mean(a, 1) + mean(a, "all");
    B = b - mean(b, 2) - mean(b, 1) + mean(b, "all");
    observed = mean(A.*B, "all");
    count = 0;
    for k = 1:nPerm
        idx = randperm(rs, n);
        count = count + (mean(A.*B(idx, idx), "all") >= observed);
    end
    pValue = (1 + count)/(nPerm + 1);
end

function q = bhfdr(p)
    p = p(:);
    n = numel(p);
    [ps, ix] = sort(p);
    adj = ps .* n ./ (1:n)';
    adj = flipud(cummin(flipud(adj)));
    q = ones(n, 1);
    q(ix) = min(adj, 1);
end
