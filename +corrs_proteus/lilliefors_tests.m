function lillieforsTests = lilliefors_tests(results)
%lilliefors_tests  Lilliefors normality test of every metric in each scenario (Lilliefors 1967, JASA 62:399).
    projectRoot = fileparts(fileparts(mfilename("fullpath")));
    conditionNames = fieldnames(results);

    alphaLevel = 0.05;
    minSample = 5;

    warnState = warning();
    warning("off", "stats:lillietest:OutOfRangePLow");
    warning("off", "stats:lillietest:OutOfRangePHigh");
    cleanup = onCleanup(@() warning(warnState));

    for field = 1:numel(conditionNames)
        conditionName = conditionNames{field};
        nMetrics = size(results.(conditionName), 2);
        outcome = nan(4, nMetrics);

        for column = 1:nMetrics
            x = results.(conditionName)(:, column);
            x = x(isfinite(x));
            if numel(x) >= minSample && std(x) > 0
                [reject, pValue, statistic, criticalValue] = lillietest(x, Alpha=alphaLevel);
                outcome(:, column) = [double(reject); pValue; statistic; criticalValue];
            end
        end

        lillieforsTests.(conditionName) = outcome;
    end

    savePath = fullfile(projectRoot, "data", 'derived', "lillieforsTests.mat");
    save(savePath, "lillieforsTests");
end
