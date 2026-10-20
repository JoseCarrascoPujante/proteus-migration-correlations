function c = coupling_theta_v(x, y, dt)
%coupling_theta_v  Speed-direction cross-correlation C_theta-v of a trajectory sampled every dt s
%   (Shaebani et al. 2022, Biophys J 121:4099, doi:10.1016/j.bpj.2022.09.033), headings measured from
%   their circular mean (Mardia & Jupp 2000, Directional Statistics, doi:10.1002/9780470316979);
%   zero-length steps, which have no heading, are left out.
    arguments
        x (:,1) double {mustBeReal}
        y (:,1) double {mustBeReal}
        dt (1,1) double {mustBePositive}
    end

    step = round(dt/0.5);   % positions are recorded every 0.5 s
    dx = diff(x(1:step:end));
    dy = diff(y(1:step:end));
    speed = hypot(dx, dy);
    % Zero-length steps have no heading (atan2(0, 0) = 0 is an artefact), so they are left out of
    % both series rather than paired with a spurious eastward heading.
    moving = speed > 0;
    speed = speed(moving);
    heading = atan2(dy(moving), dx(moving));

    headingMean = atan2(mean(sin(heading)), mean(cos(heading)));
    headingDev = mod(heading - headingMean + pi, 2*pi) - pi;
    speedDev = speed - mean(speed);

    c = mean(headingDev.*speedDev)/(std(headingDev, 1)*std(speedDev, 1));
end
