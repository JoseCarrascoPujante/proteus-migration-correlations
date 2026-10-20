function allConditionResults = fit_models(data)
%fit_models  The four two-parameter candidate forms fitted by least squares to every coupled pair in both
%   directions, with AIC and AICc.
    allConditionResults = struct();

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    metric_names = corrs_proteus.metricNames();
    varNames = metric_names;

    models = {'poly1','exp1','power1','log'};

    dep = readtable(fullfile(projectRoot, 'data', 'derived', 'pairwiseDependence.csv'), ...
        'TextType', 'string');
    isTrue = @(v) ismember(lower(string(v)), ["1", "true"]);
    keep = isTrue(dep.Dependent) & ~isTrue(dep.Definitional);
    pairNames = split(dep.Pair(keep), " vs ");
    pairNames = reshape(pairNames, [], 2);
    [~, pairs] = ismember(pairNames, string(metric_names));
    assert(all(pairs(:) > 0), 'fit_models:pairNames', 'Unknown metric name in pairwiseDependence.csv.');
    numComparisons = size(pairs, 1);

    excelFile = fullfile(projectRoot, 'data', 'derived', 'GoodnessOfFitResults.xlsx');
    if exist(excelFile, 'file'), delete(excelFile); end
    conditions = fieldnames(data);

    pooled = cellfun(@(c) data.(c), conditions, 'UniformOutput', false);
    xRef = median(vertcat(pooled{:}), 1, 'omitnan');

    for ii = 1:length(conditions)
        condition = conditions{ii};

        dataMatrix = data.(condition);

        loopResults = cell(numComparisons, 1);

        xVarSliced = cell(numComparisons, 1);
        yVarSliced = cell(numComparisons, 1);
        for jIdx = 1:numComparisons
            xVarSliced{jIdx} = varNames{pairs(jIdx, 1)};
            yVarSliced{jIdx} = varNames{pairs(jIdx, 2)};
        end

        parfor i = 1:numComparisons
            xVarIdx = pairs(i, 1);
            refs = xRef(pairs(i, :));
            yVarIdx = pairs(i, 2);

            xVar = xVarSliced{i};
            yVar = yVarSliced{i};

            x = dataMatrix(:, xVarIdx);
            y = dataMatrix(:, yVarIdx);

            gofTable = table();

            for k = 1:length(models)
                pearsonRho = corr(x, y, 'Type','Pearson', 'Rows', 'pairwise');
                spearmanRho = corr(x, y, 'Type','Spearman', 'Rows', 'pairwise');

                [gof_xy, gof_yx] = fitBothDirections(x, y, models{k}, refs);

                newRow = table(...
                    models(k), {sprintf('%s vs %s', xVar, yVar)},...
                    gof_xy.rsquare, gof_yx.rsquare, pearsonRho, ...
                    spearmanRho, gof_xy.sse, gof_yx.sse, ...
                    gof_xy.adjrsquare, gof_yx.adjrsquare, gof_xy.rmse, gof_yx.rmse, gof_xy.dfe, gof_yx.dfe,...
                    gof_xy.aic, gof_xy.aicc, gof_yx.aic, gof_yx.aicc,...
                    gof_xy.a, gof_xy.b, gof_yx.a, gof_yx.b, ...
                    gof_xy.ciA_lo, gof_xy.ciA_hi, gof_xy.ciB_lo, gof_xy.ciB_hi, ...
                    gof_yx.ciA_lo, gof_yx.ciA_hi, gof_yx.ciB_lo, gof_yx.ciB_hi, ...
                    'VariableNames', {...
                    'Model', 'Compared_Variables', ...
                    'R_squared_AxB', 'R_squared_BxA', 'Pearson_Rho', ...
                    'Spearman_Rho', 'SSE_AxB', 'SSE_BxA',...
                    'Adj_RSquared_AxB', 'Adj_RSquared_BxA', 'RMSE_AxB', 'RMSE_BxA','DFE_AxB','DFE_BxA',...
                    'AIC_AxB','AICc_AxB','AIC_BxA','AICc_BxA', ...
                    'coeff_a_AxB', 'coeff_b_AxB', 'coeff_a_BxA', 'coeff_b_BxA', ...
                    'coeff_a_lo_AxB', 'coeff_a_hi_AxB', 'coeff_b_lo_AxB', 'coeff_b_hi_AxB', ...
                    'coeff_a_lo_BxA', 'coeff_a_hi_BxA', 'coeff_b_lo_BxA', 'coeff_b_hi_BxA'});

                gofTable = [gofTable; newRow];
            end

            loopResults{i} = gofTable;
        end

        finalConditionTable = vertcat(loopResults{:});
        allConditionResults.(condition) = finalConditionTable;

        writetable(finalConditionTable, excelFile, 'Sheet', condition, 'UseExcel', false);
    end
    save(fullfile(projectRoot, 'data', 'derived', 'GoodnessOfFitResults.mat'), 'allConditionResults');
end

function [gof_xy, gof_yx] = fitBothDirections(x, y, model, refs)
    % The exponential is fitted about the pooled median of its predictor (refs), so that its scale
    % coefficient is the response at that median.
    shift = [0 0];
    if strcmp(model, 'exp1'), shift = refs; end

    [fit_xy, gof_xy] = performFit(x - shift(1), y, model);
    ic_xy = computeIC(gof_xy.sse, numcoeffs(fit_xy), numel(x));
    gof_xy.aic  = ic_xy.aic;
    gof_xy.aicc = ic_xy.aicc;
    [gof_xy.a, gof_xy.b, gof_xy.ciA_lo, gof_xy.ciA_hi, gof_xy.ciB_lo, gof_xy.ciB_hi] = getModelCoeffs(fit_xy);

    [fit_yx, gof_yx] = performFit(y - shift(2), x, model);
    ic_yx = computeIC(gof_yx.sse, numcoeffs(fit_yx), numel(y));
    gof_yx.aic  = ic_yx.aic;
    gof_yx.aicc = ic_yx.aicc;
    [gof_yx.a, gof_yx.b, gof_yx.ciA_lo, gof_yx.ciA_hi, gof_yx.ciB_lo, gof_yx.ciB_hi] = getModelCoeffs(fit_yx);
end

function [fitResult, gof] = performFit(x, y, model)
    oldWarn = warning('off', 'all');

    try
        opts = fitoptions(fittype(model));
        if strcmpi(opts.Method, 'LinearLeastSquares')
            [fitResult,gof] = fit(x, y, model);
        else
            [fitResult, gof] = fit(x, y, model,'DiffMaxChange', 0.05, ...
                    'MaxIter', 200000, 'Display', 'off');
        end
    catch ME
        warning(oldWarn);
        rethrow(ME);
    end

    warning(oldWarn);
end

function ic = computeIC(sse, kCoeff, n)
    % Gaussian log-likelihood from the residual sum of squares; k counts the coefficients and the error
    % variance, "the total number of parameters in the model, including the intercept and sigma^2"
    % (Burnham & Anderson 2004, doi:10.1177/0049124104268644).
    k = kCoeff + 1;
    const = n * (1 + log(2*pi) + log(sse/n));
    ic.aic  = const + 2*k;
    ic.aicc = ic.aic + 2*k*(k+1) / (n - k - 1);
end

function [valA, valB, ciA_lo, ciA_hi, ciB_lo, ciB_hi] = getModelCoeffs(fitResult)
    valA = NaN; valB = NaN;
    ciA_lo = NaN; ciA_hi = NaN; ciB_lo = NaN; ciB_hi = NaN;

    if isa(fitResult, 'cfit')
        cNames = coeffnames(fitResult);
        cVals = coeffvalues(fitResult);
        cInts = confint(fitResult, 0.95);

        for idx = 1:length(cNames)
            if ismember(cNames{idx}, {'a', 'p1'})
                valA = cVals(idx);
                ciA_lo = cInts(1, idx);
                ciA_hi = cInts(2, idx);
            elseif ismember(cNames{idx}, {'b', 'p2'})
                valB = cVals(idx);
                ciB_lo = cInts(1, idx);
                ciB_hi = cInts(2, idx);
            end
        end
    end
end
