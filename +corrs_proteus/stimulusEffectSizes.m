function T = stimulusEffectSizes()
%stimulusEffectSizes  Kruskal-Wallis H and epsilon-squared of each metric, and Cliff's
%   delta with its BCa CI and Mann-Whitney p of each contrast (Kruskal & Wallis 1952,
%   doi:10.1080/01621459.1952.10483441; Tomczak & Tomczak 2014, Trends Sport Sci 1:19; Mann & Whitney 1947,
%   doi:10.1214/aoms/1177730491). C_theta-v enters as its magnitude, abs_C_theta_v (see below).
    nBoot = 2000;     % bootstrap resamples of the BCa interval of Cliff's delta
    qLevel = 0.05;    % Benjamini-Hochberg threshold, also the false coverage rate of the intervals
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    D = load(fullfile(projectRoot, 'data', 'derived', 'results.mat'));
    R = D.results; conds = fieldnames(R);
    metricNames = string(corrs_proteus.metricNames());
    contrasts = [1 2; 1 3; 2 3];
    contrastNames = ["NoStim-vs-Galv","NoStim-vs-Chem","Galv-vs-Chem"];
    nMetric = numel(metricNames); nC = size(contrasts,1);
    nCell = size(R.(conds{1}), 1);
    rsBoot = RandStream("twister","Seed",3);

    Metric = strings(nMetric*nC,1); Contrast = strings(nMetric*nC,1);
    delta = zeros(nMetric*nC,1); dLo = zeros(nMetric*nC,1); dHi = zeros(nMetric*nC,1);
    boot = zeros(nBoot, nMetric*nC); z0 = zeros(nMetric*nC,1); acc = zeros(nMetric*nC,1);
    pMW = zeros(nMetric*nC,1);
    eps2 = zeros(nMetric*nC,1); Hstat = zeros(nMetric*nC,1); pKWcol = zeros(nMetric*nC,1);
    pKWmetric = zeros(nMetric,1);

    row = 0;
    for m = 1:nMetric
        G = zeros(nCell, numel(conds));
        for c = 1:numel(conds), G(:,c) = R.(conds{c})(:, m); end
        name = metricNames(m);
        % A mirror image of a trajectory reverses the sign of C_theta-v and leaves every other metric
        % unchanged, so the scenarios are compared on |C_theta-v|, which a reflection leaves unchanged too.
        if name == "C_theta_v"
            G = abs(G);
            name = "abs_C_theta_v";
        end
        [Hm, pKW, e2] = kwOmnibus(G);
        pKWmetric(m) = pKW;
        for k = 1:nC
            row = row + 1;
            xi = G(:, contrasts(k,1)); yj = G(:, contrasts(k,2));
            Metric(row) = name; Contrast(row) = contrastNames(k);
            delta(row) = corrs_proteus.cliffsDelta(xi, yj);
            [boot(:,row), z0(row), acc(row)] = bcaCliffBoot(xi, yj, nBoot, rsBoot);
            pMW(row) = ranksum(xi, yj);
            eps2(row) = e2; Hstat(row) = Hm; pKWcol(row) = pKW;
        end
    end

    qMW = bhfdr(pMW);
    % The intervals are drawn after the Benjamini-Hochberg selection, at the level 1 - R*qLevel/m that bounds
    % the false coverage rate of the R contrasts declared significant among the m tested by qLevel
    % (Benjamini & Yekutieli 2005, doi:10.1198/016214504000001907); the bootstrap draws are those above.
    % Here 12 of the 33 contrasts have q < 0.05, so the level is 98.2% (Materials and Methods, Cliff's delta).
    alphaCI = sum(qMW < qLevel) * qLevel / numel(qMW);
    for r = 1:numel(qMW)
        [dLo(r), dHi(r)] = bcaLimits(boot(:,r), z0(r), acc(r), alphaCI);
    end
    qKWmetric = bhfdr(pKWmetric);
    KW_qFDR = repelem(qKWmetric, nC);

    T = table(Metric, Contrast, delta, dLo, dHi, pMW, qMW, eps2, Hstat, pKWcol, KW_qFDR, ...
        'VariableNames', {'Metric','Contrast','CliffsDelta','delta_lo','delta_hi', ...
        'MannWhitney_p','MannWhitney_qFDR','KW_epsilon2','KW_H','KW_p','KW_qFDR'});
    for vn = T.Properties.VariableNames
        if isnumeric(T.(vn{1})), T.(vn{1}) = corrs_proteus.roundSigMin3(T.(vn{1})); end
    end
    out = fullfile(projectRoot, 'data', 'derived', 'stimulusEffectSizes.csv');
    writetable(T, out);
end

function [H, p, eps2] = kwOmnibus(G)
    [p, tbl] = kruskalwallis(G, [], 'off');
    H = tbl{2, 5};
    eps2 = H / (numel(G) - 1);
end

function [boot, z0, a] = bcaCliffBoot(x, y, B, rs)
    nx = numel(x); ny = numel(y);
    theta = corrs_proteus.cliffsDelta(x, y);
    boot = zeros(B,1);
    for b = 1:B
        boot(b) = corrs_proteus.cliffsDelta(x(randi(rs,nx,nx,1)), y(randi(rs,ny,ny,1)));
    end
    z0 = norminv(mean(boot < theta));
    jack = zeros(nx+ny,1); t = 0;
    for i = 1:nx, t=t+1; jack(t) = corrs_proteus.cliffsDelta(x([1:i-1,i+1:nx]), y); end
    for j = 1:ny, t=t+1; jack(t) = corrs_proteus.cliffsDelta(x, y([1:j-1,j+1:ny])); end
    dd = mean(jack) - jack;
    den = 6 * (sum(dd.^2))^1.5;
    a = sum(dd.^3) / den;
end

function [lo, hi] = bcaLimits(boot, z0, a, alpha)
    za = norminv(alpha/2); z1 = norminv(1 - alpha/2);
    a1 = normcdf(z0 + (z0 + za) / (1 - a*(z0 + za)));
    a2 = normcdf(z0 + (z0 + z1) / (1 - a*(z0 + z1)));
    lo = quantile(boot, a1); hi = quantile(boot, a2);
end

function q = bhfdr(p)
    p = p(:); n = numel(p); [ps, ix] = sort(p);
    adj = ps .* n ./ (1:n)'; adj = flipud(cummin(flipud(adj)));
    q = ones(n,1); q(ix) = min(adj, 1);
end
