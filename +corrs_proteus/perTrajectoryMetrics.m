function perTrajectoryMetrics()
%perTrajectoryMetrics  The value of each metric in every cell, one table per metric.
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    dataDir = fullfile(projectRoot, 'data', 'derived');
    D = load(fullfile(dataDir, 'results.mat'));
    R = D.results; conds = string(fieldnames(R))';

    metricNames = string(corrs_proteus.metricNames());
    condDisp = mapCondNames(conds);
    nCell = size(R.(conds(1)), 1);

    for m = 1:numel(metricNames)
        vals = zeros(nCell, numel(conds));
        for c = 1:numel(conds)
            vals(:, c) = R.(conds(c))(:, m);
        end
        if metricNames(m) == "RMSF_correlation_time"
            vals = vals / 120;
        end
        vals = round(vals, 3);
        T = array2table([(1:nCell)', vals], 'VariableNames', ["Cell", condDisp]);
        out = fullfile(dataDir, sprintf('perTrajectory_%s.csv', metricNames(m)));
        writetable(T, out);
    end
end

function out = mapCondNames(conds)
    keys = ["NoStimuli","Galvanotaxis","Chemotaxis"];
    vals = ["No_stimuli","Galvanotaxis","Chemotaxis"];
    out = string(conds);
    [tf, loc] = ismember(out, keys);
    out(tf) = vals(loc(tf));
end
