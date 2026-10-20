function T = pairwiseDependence()
%pairwiseDependence  Dependence of the 55 metric pairs.
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    D = load(fullfile(projectRoot, 'data', 'derived', 'results.mat'));
    T = corrs_proteus.dependence_rule(D.results);
    % Rows by Rank (decreasing median R*_n), so that every table built from this file shares the order.
    T = sortrows(T, 'Rank');
    out = fullfile(projectRoot, 'data', 'derived', 'pairwiseDependence.csv');
    writetable(T, out);
end
