function summary = perContestDecisionSummary(bordaResult, modelFitResults)
%perContestDecisionSummary  AICc ranks and fitted parameters, one row per contest.
    arguments
        bordaResult (1,1) struct
        modelFitResults (1,1) struct
    end
    focusModel = "exp1";   % the candidate ranked first most often by AICc

    contestKeys = {'Condition', 'Pair', 'Direction'};
    L = bordaResult.longTable;
    fm = char(focusModel);

    conditions = fieldnames(modelFitResults);
    flatData = table();

    for c = 1:numel(conditions)
        condName = conditions{c};
        tbl = modelFitResults.(condName);


        for direction = ["AxB", "BxA"]
            dirStr = char(direction);

            subTbl = table();
            subTbl.Condition = repmat(string(condName), height(tbl), 1);
            subTbl.Pair      = string(tbl.Compared_Variables);
            subTbl.Direction = repmat(direction, height(tbl), 1);
            subTbl.Model     = string(tbl.Model);

            subTbl.Rsquare = tbl.(['R_squared_' dirStr]);
            subTbl.a       = tbl.(['coeff_a_' dirStr]);
            subTbl.b       = tbl.(['coeff_b_' dirStr]);
            subTbl.a_lo    = tbl.(['coeff_a_lo_' dirStr]);
            subTbl.a_hi    = tbl.(['coeff_a_hi_' dirStr]);
            subTbl.b_lo    = tbl.(['coeff_b_lo_' dirStr]);
            subTbl.b_hi    = tbl.(['coeff_b_hi_' dirStr]);

            flatData = [flatData; subTbl];
        end
    end

    winners = L(L.rank == 1, [contestKeys, {'Model','AICc','akaikeWeight'}]);
    winners.Properties.VariableNames(end-2:end) = ...
        {'WinnerModel','WinnerAICc','WinnerAkaikeWeight'};

    F = L(L.Model == focusModel, ...
        [contestKeys, {'AICc','deltaAICc','akaikeWeight'}]);
    F.Properties.VariableNames(end-2:end) = ...
        {'FocusAICc','FocusDeltaAICc','FocusAkaikeWeight'};

    summary = outerjoin(winners, F, 'Keys', contestKeys, ...
        'MergeKeys', true, 'Type', 'left');

    winR2 = flatData(:, [contestKeys, {'Model', 'Rsquare'}]);
    winR2.Properties.VariableNames{'Model'} = 'WinnerModel';
    winR2.Properties.VariableNames{'Rsquare'} = 'WinnerR2';
    summary = outerjoin(summary, winR2, ...
        'Keys', [contestKeys, {'WinnerModel'}], 'MergeKeys', true, 'Type', 'left');

    focusParams = flatData(flatData.Model == focusModel, ...
        [contestKeys, {'Rsquare', 'a','a_lo','a_hi', 'b','b_lo','b_hi'}]);

    focusParams = renamevars(focusParams, ...
        {'Rsquare','a','a_lo','a_hi','b','b_lo','b_hi'}, ...
        {'FocusR2','FocusA','FocusA_lo','FocusA_hi','FocusB','FocusB_lo','FocusB_hi'});

    summary = outerjoin(summary, focusParams, 'Keys', contestKeys, ...
        'MergeKeys', true, 'Type', 'left');

    summary = sortrows(summary, {'Condition','Pair','Direction'});
    summary = renamevars(summary, ...
        {'FocusR2','FocusAICc','FocusDeltaAICc','FocusAkaikeWeight', ...
         'FocusA','FocusA_lo','FocusA_hi', ...
         'FocusB','FocusB_lo','FocusB_hi'}, ...
        compose("%s_%s", fm, ...
            ["R2","AICc","DeltaAICc","AkaikeWeight", ...
             "a","a_lo","a_hi","b","b_lo","b_hi"]));

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    writetable(summary, fullfile(projectRoot, 'data', 'derived', sprintf('perContestDecision_%s.csv', fm)));
end
