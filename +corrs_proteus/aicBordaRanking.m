function results = aicBordaRanking()
%aicBordaRanking  Borda count of the AICc rankings of the candidate forms over all contests.
%   AICc: Akaike 1974, doi:10.1109/TAC.1974.1100705; Sugiura 1978, doi:10.1080/03610927808827599;
%   Hurvich & Tsai 1989, doi:10.1093/biomet/76.2.297. Akaike weights: Burnham & Anderson 2004,
%   doi:10.1177/0049124104268644. Borda count: Brams & Fishburn 2002, doi:10.1016/S1574-0110(02)80008-X;
%   Blesch et al. 2022, doi:10.1038/s41562-022-01430-7.
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    F = load(fullfile(projectRoot, 'data', 'derived', 'GoodnessOfFitResults.mat'), 'allConditionResults');
    conds = fieldnames(F.allConditionResults);
    allRows = table();
    for i = 1:numel(conds)
        c = F.allConditionResults.(conds{i});
        c.Condition = repmat(string(conds{i}), height(c), 1);
        allRows = [allRows; c]; %#ok<AGROW>
    end

    % Both fit directions of every pair and scenario are contests.
    longTable = [ ...
        table(repmat("AxB", height(allRows), 1), string(allRows.Condition), ...
            string(allRows.Compared_Variables), string(allRows.Model), allRows.AICc_AxB, ...
            'VariableNames', {'Direction','Condition','Pair','Model','AICc'}); ...
        table(repmat("BxA", height(allRows), 1), string(allRows.Condition), ...
            string(allRows.Compared_Variables), string(allRows.Model), allRows.AICc_BxA, ...
            'VariableNames', {'Direction','Condition','Pair','Model','AICc'})];

    longTable.ContestID = longTable.Condition + "|" + longTable.Pair + "|" + longTable.Direction;

    G = findgroups(longTable.ContestID);
    longTable.deltaAICc    = nan(height(longTable),1);
    longTable.akaikeWeight = nan(height(longTable),1);
    longTable.rank         = nan(height(longTable),1);
    longTable.bordaPts     = nan(height(longTable),1);

    for g = 1:max(G)
        idx  = find(G == g);
        aicc = longTable.AICc(idx);
        d    = aicc - min(aicc);
        w    = exp(-d/2) / sum(exp(-d/2));
        [~, ord] = sort(aicc, 'ascend');
        r        = zeros(size(aicc));
        r(ord)   = 1:numel(aicc);
        n        = numel(aicc);
        b        = n - r;

        longTable.deltaAICc(idx)    = d;
        longTable.akaikeWeight(idx) = w;
        longTable.rank(idx)         = r;
        longTable.bordaPts(idx)     = b;
    end

    [g, modelsG]   = findgroups(longTable.Model);
    totalBorda     = splitapply(@sum,   longTable.bordaPts,     g);
    contestsEntered   = splitapply(@numel, longTable.bordaPts,     g);
    meanAICc       = splitapply(@mean,  longTable.AICc,         g);
    meanDelta      = splitapply(@mean,  longTable.deltaAICc,    g);
    meanWeight     = splitapply(@mean,  longTable.akaikeWeight, g);
    meanRank       = splitapply(@mean,  longTable.rank,         g);
    firstPlaces    = splitapply(@(r) sum(r == 1), longTable.rank, g);
    nContestsTotal    = numel(unique(longTable.ContestID));

    globalTbl = table( ...
        modelsG, totalBorda, firstPlaces, contestsEntered, ...
        repmat(nContestsTotal, numel(modelsG), 1), ...
        totalBorda ./ contestsEntered, ...
        meanRank, meanDelta, meanAICc, meanWeight, ...
        'VariableNames', {'Model','TotalBorda','FirstPlaces','ContestsEntered','TotalContests', ...
                          'MeanBordaPerContest','MeanRank','MeanDeltaAICc', ...
                          'MeanAICc','MeanAkaikeWeight'});
    globalTbl = sortrows(globalTbl, 'TotalBorda', 'descend');

    [g, modelsC, condsC] = findgroups(longTable.Model, longTable.Condition);
    totalBordaC = splitapply(@sum,  longTable.bordaPts,     g);
    weightC     = splitapply(@mean, longTable.akaikeWeight, g);
    perCondTbl  = table(modelsC, condsC, totalBordaC, weightC, ...
        'VariableNames', {'Model','Condition','TotalBorda','MeanAkaikeWeight'});
    perCondTbl  = sortrows(perCondTbl, {'Condition','TotalBorda'}, {'ascend','descend'});

    [g, modelsP, pairsP] = findgroups(longTable.Model, longTable.Pair);
    totalBordaP = splitapply(@sum,  longTable.bordaPts,     g);
    weightP     = splitapply(@mean, longTable.akaikeWeight, g);
    perPairTbl  = table(modelsP, pairsP, totalBordaP, weightP, ...
        'VariableNames', {'Model','Pair','TotalBorda','MeanAkaikeWeight'});
    perPairTbl  = sortrows(perPairTbl, {'Pair','TotalBorda'}, {'ascend','descend'});

    dataDir = fullfile(projectRoot, 'data', 'derived');
    writetable(globalTbl,  fullfile(dataDir, 'aicBorda_global.csv'));
    writetable(perCondTbl, fullfile(dataDir, 'aicBorda_perCondition.csv'));
    writetable(perPairTbl, fullfile(dataDir, 'aicBorda_perPair.csv'));
    writetable(longTable,  fullfile(dataDir, 'aicBorda_longTable.csv'));

    results = struct( ...
        'global',       globalTbl, ...
        'perCondition', perCondTbl, ...
        'perPair',      perPairTbl, ...
        'longTable',    longTable);
end
