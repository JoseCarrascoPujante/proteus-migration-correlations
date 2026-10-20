function aiccBordaRankingTable()
%aiccBordaRankingTable  AICc Borda ranking of the candidate forms.
    dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'derived');
    b = readtable(fullfile(dataDir, "aicBorda_global.csv"), TextType = "string");
    [~, ord] = sort(b.MeanBordaPerContest, "descend");
    rnk = zeros(height(b), 1); rnk(ord) = 1:height(b); b.Rank = rnk;
    s2 = b(:, {'Rank','Model','MeanBordaPerContest','ContestsEntered', ...
        'TotalBorda','FirstPlaces','MeanDeltaAICc','MeanAkaikeWeight','MeanRank'});
    s2 = sortrows(s2, 'Rank');
    for vn = s2.Properties.VariableNames
        if isnumeric(s2.(vn{1})), s2.(vn{1}) = corrs_proteus.roundSigMin3(s2.(vn{1})); end
    end
    out = fullfile(dataDir, 'aiccBordaRanking.csv');
    writetable(s2, out);
end
