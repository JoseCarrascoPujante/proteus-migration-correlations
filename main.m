% Reproduces every result and table of the article from the trajectories in data/coordinates.mat.
% Run as main from this folder (see README.md).

clearvars; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
if isempty(projectRoot), projectRoot = pwd; end
derivedDir = fullfile(projectRoot, 'data', 'derived');

addpath(projectRoot);
if ~isfolder(derivedDir), mkdir(derivedDir); end

%% 1. Trajectories: x and y (mm) of every cell each 0.5 s, from its initial position
load(fullfile(projectRoot, 'data', 'coordinates.mat'), 'coordinates');

%% 2. Metrics of every cell
results = corrs_proteus.metric_calculation(coordinates, derivedDir);

%% 3. Lilliefors normality tests
lillieforsResults = corrs_proteus.lilliefors_tests(results);

%% 4. Dependence of the 55 metric pairs
corrs_proteus.pairwiseDependence();

%% 5. Stimulus effect sizes
corrs_proteus.stimulusEffectSizes();

%% 6. Between-scenario differences in the dependence of every metric pair
corrs_proteus.dependenceContrasts();

%% 7. Candidate forms fitted to every dependent pair
modelFitResults = corrs_proteus.fit_models(results);

%% 8. AICc ranking and Borda count
resOLS = corrs_proteus.aicBordaRanking();

%% 9. AICc ranks and parameters of every contest
summaryOLS = corrs_proteus.perContestDecisionSummary(resOLS, modelFitResults);

%% 10. Strength and direction of the six couplings
couplings = corrs_proteus.couplingStrength();

%% 11. AICc Borda ranking of the candidate forms
corrs_proteus.aiccBordaRankingTable();

%% 12. Value of each metric in every cell
corrs_proteus.perTrajectoryMetrics();

%% 13. Departure of DFA gamma from the self-similar value beta/2 + 1
iGamma = find(strcmp(corrs_proteus.metricNames(), 'DFA_gamma'));
iBeta = find(strcmp(corrs_proteus.metricNames(), 'MSD_beta'));
for cond = string(fieldnames(results))'
    M = results.(cond);
    departure = abs(M(:, iGamma) - (M(:, iBeta) / 2 + 1));
    fprintf('%-13s median |gamma - (beta/2 + 1)| = %.3f   IQR(gamma) = %.3f\n', ...
        cond, median(departure), iqr(M(:, iGamma)));
end

%% 14. Robustness of the every-frame metrics to the zero-length steps
% Manual tracking records the clicked pixel, so a cell that moves less than one pixel between frames leaves a
% zero-length step (Gorelik & Gautreau 2014, doi:10.1038/nprot.2014.131). The path length behind DR and AS,
% and the RMSF exponent, are recomputed from the positions at the resolved interval, where steps are almost
% never zero, and compared with their values from every frame.
kRes = round(corrs_proteus.resolvedInterval(coordinates) / 0.5);
names = string(corrs_proteus.metricNames());
col = @(n) find(names == n);
dc = @(a, b) corrs_proteus.distance_correlation(a, b);
for cond = string(fieldnames(coordinates))'
    X = coordinates.(cond).x; Y = coordinates.(cond).y; M = results.(cond);
    d = hypot(diff(X), diff(Y));
    dRes = hypot(diff(X(1:kRes:end, :)), diff(Y(1:kRes:end, :)));
    drRes = M(:, col("IR")) ./ sum(dRes, 1)';
    asRes = sum(dRes, 1)' * 1000 / 1800;
    alphaRes = arrayfun(@(j) corrs_proteus.rmsf(dRes(:, j), kRes * 0.5), (1:size(X, 2))');
    fprintf(['%-13s zero steps %.1f%%; L(0.5 s)/L(%.1f s) median %.3f; ' ...
        'Spearman with the every-frame values: DR %.3f, AS %.3f, RMSF alpha %.3f\n'], cond, ...
        100 * mean(d(:) == 0), kRes * 0.5, median(sum(d, 1)' ./ sum(dRes, 1)'), ...
        corr(M(:, col("DR")), drRes, 'Type', 'Spearman'), corr(M(:, col("AS")), asRes, 'Type', 'Spearman'), ...
        corr(M(:, col("RMSF_alpha")), alphaRes, 'Type', 'Spearman'));
    fprintf(['%-13s R*_n at the resolved interval (every frame): DFA gamma-DR %.3f (%.3f), ' ...
        'MSD beta-DR %.3f (%.3f), RMSF alpha-AC_v(1) %.3f (%.3f)\n'], cond, ...
        dc(M(:, col("DFA_gamma")), drRes), dc(M(:, col("DFA_gamma")), M(:, col("DR"))), ...
        dc(M(:, col("MSD_beta")), drRes), dc(M(:, col("MSD_beta")), M(:, col("DR"))), ...
        dc(alphaRes, M(:, col("AC_v1"))), dc(M(:, col("RMSF_alpha")), M(:, col("AC_v1"))));
end

%% 15. Side of the net displacement
% The cathode lies to the right (+x) in the galvanotaxis set-up, and the nFMLP was added on the left (-x) in the
% chemotaxis set-up; each cell is counted on the side toward which its final position lies from its initial one.
for cond = string(fieldnames(coordinates))'
    netX = coordinates.(cond).x(end, :) - coordinates.(cond).x(1, :);
    fprintf('%-13s net displacement to the right (+x): %d cells; to the left (-x): %d cells; of %d\n', ...
        cond, nnz(netX > 0), nnz(netX < 0), numel(netX));
end

%% 16. Speed memory against temporally independent speeds
% For n independent values the lag-one sample autocorrelation is approximately normal with mean 0 and variance
% 1/n (Brockwell & Davis 2016, Section 1.6, doi:10.1007/978-3-319-29854-2), so its two-sided 95% bound is the
% 0.975 normal quantile over sqrt(n), n being the number of speeds of each track at the resolved interval.
conds = fieldnames(coordinates);
nSpeeds = numel(1:kRes:size(coordinates.(conds{1}).x, 1)) - 1;
acv = cell2mat(cellfun(@(c) results.(c)(:, col("AC_v1")), conds, 'UniformOutput', false));
fprintf('AC_v(1) from %.3f to %.3f over the %d cells; 95%% bound for independent speeds over %d steps: +/-%.3f\n', ...
    min(acv), max(acv), numel(acv), nSpeeds, norminv(0.975) / sqrt(nSpeeds));

%% 17. Sampling-interval sweep
sweep = corrs_proteus.samplingIntervalSweep(coordinates, results);
% The values the article quotes: R*_n of AC_v(1) with RMSF alpha from the resolved interval to 60 s, the R*_n of
% t_theta at the intervals of 15 s and longer, and that of t_theta with AS at any interval.
kLong = round(max(sweep.dt_s) / 0.5);          % positions are recorded every 0.5 s
fprintf('Steps per track: %d at %.1f s, %d at %.0f s\n', nSpeeds, kRes * 0.5, ...
    numel(1:kLong:size(coordinates.(conds{1}).x, 1)) - 1, kLong * 0.5);
for cond = string(conds)'
    S = sweep(sweep.Scenario == cond, :);
    fromRes = S.dt_s >= kRes * 0.5;
    long = S.dt_s >= 15;
    fprintf(['%-13s R*_n RMSF alpha-AC_v(1) from %.1f s: min %.3f; from 15 s: t_theta-MSD beta %.3f to %.3f, ' ...
        't_theta-AS %.3f to %.3f; t_theta-AS at any interval: max %.3f\n'], cond, kRes * 0.5, ...
        min(S.bcDcor_RMSF_alpha_AC_v1(fromRes)), min(S.bcDcor_t_theta_MSD_beta(long)), ...
        max(S.bcDcor_t_theta_MSD_beta(long)), min(S.bcDcor_t_theta_AS(long)), max(S.bcDcor_t_theta_AS(long)), ...
        max(S.bcDcor_t_theta_AS));
end

%% 18. Principal component analysis of the metrics
[~, ~, prTable] = corrs_proteus.principalComponents(results);
fprintf('Participation ratio: %s\n', strjoin(compose("%s %.2f", prTable.Scope, prTable.ParticipationRatio), ', '));

%% 19. Dendrogram of the metrics by their dependence
corrs_proteus.dependenceDendrogram(results);

%% 20. Two-lines U-shape test (R)
% Rscript is looked for on the system path and then in the default Windows installation folder.
rscript = "Rscript";
[notOnPath, ~] = system('Rscript --version');
if notOnPath
    found = dir('C:\Program Files\R\R-*\bin\Rscript.exe');
    assert(~isempty(found), 'main:noR', 'R not found; the two-lines test needs it (see README.md).');
    rscript = string(fullfile(found(end).folder, found(end).name));
end
status = system(sprintf('"%s" "%s" "%s"', rscript, fullfile(projectRoot, 'R', 'twolines_check.R'), projectRoot));
assert(status == 0, 'main:twolines', 'The two-lines test failed; see the R output above.');
