function T = samplingIntervalSweep(coordinates, results)
%samplingIntervalSweep  Fig S4: the four metrics computed at the resolved interval, and the dependences
%   involving them on which the paper draws conclusions, at every sampling interval from 0.5 to 60 s.
%   Only ApEn, AC_v(1), t_theta and C_theta-v have an interval to choose, so only they are swept. For every
%   interval Dt and scenario, T holds the 25th, 50th and 75th percentiles over the cells of each of the four
%   metrics, computed as in metric_calculation (cells whose mean cosine of the turning angle is not positive
%   have no t_theta at that interval and are left out), and the bias-corrected squared distance correlation
%   R*_n of AC_v(1) with RMSF alpha and of t_theta with MSD beta and with AS. RMSF alpha, MSD beta and AS use
%   every frame and keep their values in results, so only AC_v(1) and t_theta change with Dt.
%   Writes data/derived/samplingIntervalSweep.csv.
    arguments
        coordinates (1, 1) struct
        results (1, 1) struct
    end
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    names = string(corrs_proteus.metricNames());
    iAlpha = find(names == "RMSF_alpha");
    partner = [find(names == "MSD_beta"), find(names == "AS")];
    steps = 1:120;                                  % positions are recorded every 0.5 s: 0.5 to 60 s
    dts = steps * 0.5;
    conds = fieldnames(coordinates);
    swept = ["ApEn", "AC_v1", "t_theta", "C_theta_v"];
    nRow = numel(conds) * numel(steps);
    [Scenario, dt_s] = deal(strings(nRow, 1), zeros(nRow, 1));
    Q = zeros(nRow, 3, numel(swept));               % row x [q25 median q75] x swept metric
    R = zeros(nRow, 1 + numel(partner));            % row x dependence
    row = 0;
    for c = 1:numel(conds)
        X = coordinates.(conds{c}).x;
        Y = coordinates.(conds{c}).y;
        M = results.(conds{c});
        nCell = size(X, 2);
        for s = 1:numel(steps)
            row = row + 1;
            Scenario(row) = conds{c};
            dt_s(row) = dts(s);
            v = nan(nCell, numel(swept));
            for j = 1:nCell
                us = hypot(diff(X(1:steps(s):end, j)), diff(Y(1:steps(s):end, j)));
                v(j, 1) = corrs_proteus.ApEn(2, 0.2 * std(us), us);
                v(j, 2) = corrs_proteus.speed_memory(X(:, j), Y(:, j), dts(s));
                v(j, 3) = corrs_proteus.direction_memory(X(:, j), Y(:, j), dts(s));
                v(j, 4) = corrs_proteus.coupling_theta_v(X(:, j), Y(:, j), dts(s));
            end
            for k = 1:numel(swept)
                ok = ~isnan(v(:, k));
                Q(row, :, k) = quantile(v(ok, k), [0.25 0.5 0.75]);
            end
            R(row, 1) = corrs_proteus.distance_correlation(M(:, iAlpha), v(:, 2));
            ok = ~isnan(v(:, 3));
            for k = 1:numel(partner)
                R(row, 1 + k) = corrs_proteus.distance_correlation(M(ok, partner(k)), v(ok, 3));
            end
        end
    end
    T = table(Scenario, dt_s);
    for k = 1:numel(swept)
        T.(swept(k) + "_q25") = Q(:, 1, k);
        T.(swept(k) + "_median") = Q(:, 2, k);
        T.(swept(k) + "_q75") = Q(:, 3, k);
    end
    T.bcDcor_RMSF_alpha_AC_v1 = R(:, 1);
    T.bcDcor_t_theta_MSD_beta = R(:, 2);
    T.bcDcor_t_theta_AS = R(:, 3);
    writetable(T, fullfile(projectRoot, 'data', 'derived', 'samplingIntervalSweep.csv'));
end
