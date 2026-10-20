function [tTheta, cosAlpha] = direction_memory(x, y, dt)
%direction_memory  Directional persistence time t_theta = -dt/ln<cos alpha> of a trajectory sampled every dt s
%   (Shaebani et al. 2022, Biophys J 121:4099, doi:10.1016/j.bpj.2022.09.033); zero-length steps,
%   which have no heading, are left out.
    arguments
        x (:,1) double {mustBeReal}
        y (:,1) double {mustBeReal}
        dt (1,1) double {mustBePositive}
    end

    step = round(dt/0.5);   % positions are recorded every 0.5 s
    dx = diff(x(1:step:end));
    dy = diff(y(1:step:end));
    % A zero-length step has no heading: atan2(0, 0) = 0 would
    % insert a spurious eastward heading and two spurious turns. Such steps are left out, so each
    % turning angle is taken between consecutive steps that moved.
    moving = dx ~= 0 | dy ~= 0;
    heading = atan2(dy(moving), dx(moving));
    cosAlpha = mean(cos(diff(heading)));

    if cosAlpha > 0
        tTheta = -dt/log(cosAlpha);
    else
        tTheta = NaN;
    end
end
