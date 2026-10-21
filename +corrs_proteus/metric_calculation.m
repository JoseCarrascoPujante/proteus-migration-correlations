function results = metric_calculation(coordinates, destination_folder)
%metric_calculation  The eleven metrics of every trajectory.
    condition = fieldnames(coordinates) ;

    stat_names = corrs_proteus.metricNames() ;

    dtMemory = corrs_proteus.resolvedInterval(coordinates) ;
    kMemory = round(dtMemory/0.5) ;
    fprintf('Resolved interval: %.1f s (%d steps per track)\n', dtMemory, ...
        floor((size(coordinates.(condition{1}).x, 1) - 1)/kMemory)) ;

    for i = 1:length(condition)

        N = size(coordinates.(condition{i}).x, 2);
        res_cond = zeros(N, length(stat_names)) ;

        parfor j = 1:N

            row = zeros(1, numel(stat_names)) ;

            xj = coordinates.(condition{i}).x(:,j);
            yj = coordinates.(condition{i}).y(:,j);
            u = hypot(diff(xj), diff(yj));

            [row(1), row(2)] = corrs_proteus.rmsf(u, 0.5) ;

            row(3) = corrs_proteus.DFA(hypot(xj, yj)) ;

            row(4) = corrs_proteus.msd(xj, yj) ;

            us = hypot(diff(xj(1:kMemory:end)), diff(yj(1:kMemory:end)));
            % m = 2, the run length at which ApEn can be estimated from relatively few points, and r = 0.2 SD of
            % the series (Pincus 1991, doi:10.1073/pnas.88.6.2297).
            row(5) = corrs_proteus.ApEn(2, 0.2*std(us), us) ;

            row(6) = ...
                norm([coordinates.(condition{i}).x(end,j) coordinates.(condition{i}).y(end,j)]...
                - [coordinates.(condition{i}).x(1,j) coordinates.(condition{i}).y(1,j)]);

            d = hypot(diff(coordinates.(condition{i}).x(:,j)), diff(coordinates.(condition{i}).y(:,j))) ;
            distTrav = sum(d);
            row(7) = row(6)/distTrav;

            % AS = path length over the duration of the recording, 3,600 frames of 0.5 s,
            % the T of the identity IR = DR*AS*T; x 1000 for mm to um.
            row(8) = distTrav*1000/(numel(xj)*0.5);

            row(9)  = corrs_proteus.coupling_theta_v(xj, yj, dtMemory);
            row(10) = corrs_proteus.direction_memory(xj, yj, dtMemory);
            row(11) = corrs_proteus.speed_memory(xj, yj, dtMemory);

            res_cond(j,:) = row ;

        end

        results.(condition{i}) = res_cond ;
    end

    assert(all(structfun(@(M) all(isfinite(M), 'all'), results)), 'metric_calculation:missing', ...
        'A metric is NaN or Inf for some cell; resolve it before any analysis.');

    save(fullfile(destination_folder, 'results.mat'), 'results') ;
end
