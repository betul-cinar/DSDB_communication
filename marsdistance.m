%==========================================================================
% worst_case_range_mission.m
% Compute global max/min SC-Earth distance for a Mars orbiter at 300 km altitude
% over the mission period (01 Nov 2031 to 31 Dec 2033), and visualize geometry.
%==========================================================================

clear; clc; close all

%% 1) Orbit parameters
R_mars  = 3396.2;           % Mars radius [km]
alt     = 300;              % Orbiter altitude above surface [km]
r_orb   = R_mars + alt;     % Orbital radius [km]

incl  = 0;   RAAN  = 0;     argp  = 0;   f0    = 0;

%% 2) Mission epoch vector
epochs = (datetime(2031,11,1):days(1):datetime(2033,12,31))';
Ne = numel(epochs);

%% 3) Pre-sample orbit (true anomaly)
Ntheta = 360;
theta  = linspace(0,2*pi,Ntheta);
r_pf   = [r_orb*cos(theta); r_orb*sin(theta); zeros(1,Ntheta)];
% rotation to Mars-ECI
Rz     = @(phi)[cosd(phi) -sind(phi) 0; sind(phi) cosd(phi) 0; 0 0 1];
rot    = Rz(RAAN) * [1 0 0; 0 cosd(incl) -sind(incl); 0 sind(incl) cosd(incl)] * Rz(argp);

%% 4) Find global extremes
globalMax = -Inf; globalMin = Inf;
maxInfo = struct('range',0,'epoch',epochs(1),'anom',0);
minInfo = struct('range',0,'epoch',epochs(1),'anom',0);
for k = 1:Ne
    JD       = juliandate(epochs(k));
    rSun2E   = planetEphemeris(JD,'Sun','Earth')';
    rSun2M   = planetEphemeris(JD,'Sun','Mars')';
    rMars2E  = rSun2E - rSun2M;
    rOrbitECI= rot * r_pf;
    % spacecraft inertial
    rSun2SC  = rSun2M + rOrbitECI;
    rSC2E    = rSun2E - rSun2SC;
    ranges   = vecnorm(rSC2E,2,1);
    [rMax,iMax] = max(ranges);
    [rMin,iMin] = min(ranges);
    if rMax > globalMax  % update max
        globalMax     = rMax;
        maxInfo.range = rMax;
        maxInfo.epoch = epochs(k);
        maxInfo.anom  = rad2deg(theta(iMax));
    end
    if rMin < globalMin  % update min
        globalMin     = rMin;
        minInfo.range = rMin;
        minInfo.epoch = epochs(k);
        minInfo.anom  = rad2deg(theta(iMin));
    end
end

%% 5) Report extremes
fprintf('Global MAX range = %.1f km  at %s, anomaly = %.1f°\n', ...
    maxInfo.range, datestr(maxInfo.epoch,'dd-mmm-yyyy'), maxInfo.anom);
fprintf('Global MIN range = %.1f km  at %s, anomaly = %.1f°\n', ...
    minInfo.range, datestr(minInfo.epoch,'dd-mmm-yyyy'), minInfo.anom);

%% 6) Recompute positions at extremes
% Max
thetaMax     = deg2rad(maxInfo.anom);
rOrbitECI_max= rot * [r_orb*cos(thetaMax); r_orb*sin(thetaMax); 0];
JDmax        = juliandate(maxInfo.epoch);
rSun2E_max   = planetEphemeris(JDmax,'Sun','Earth')';
rSun2M_max   = planetEphemeris(JDmax,'Sun','Mars')';
rMars2E_max  = rSun2E_max - rSun2M_max;
% Min
thetaMin     = deg2rad(minInfo.anom);
rOrbitECI_min= rot * [r_orb*cos(thetaMin); r_orb*sin(thetaMin); 0];
JDmin        = juliandate(minInfo.epoch);
rSun2E_min   = planetEphemeris(JDmin,'Sun','Earth')';
rSun2M_min   = planetEphemeris(JDmin,'Sun','Mars')';
rMars2E_min  = rSun2E_min - rSun2M_min;

%% 7) Plot: two subplots for clarity
figure('Name','Range Extremes Visualization','Position',[100 100 1200 500]);

%---- Left: Local Mars view ----
subplot(1,2,1);
hold on; axis equal; grid on;
% Mars
c = linspace(0,2*pi,200);
fill(R_mars*cos(c), R_mars*sin(c), [1 .6 .4], 'FaceAlpha', .5, 'EdgeColor', 'none');
plot(R_mars*cos(c), R_mars*sin(c), 'k-', 'LineWidth', 1);
% Orbit path
p = linspace(0,2*pi,400);
plot(r_orb*cos(p), r_orb*sin(p), 'r--', 'LineWidth', 1);
% Spacecraft markers filled
plot(rOrbitECI_max(1), rOrbitECI_max(2), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'SC @ MAX');
plot(rOrbitECI_min(1), rOrbitECI_min(2), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'SC @ MIN');
legend('Mars Surface', 'Mars Outline', 'Orbit Path', 'SC @ MAX', 'SC @ MIN', 'Location', 'best');
title('Mars Orbit Zoom');
xlabel('X [km]'); ylabel('Y [km]');

%---- Right: SC-Earth vectors ----
subplot(1,2,2);
hold on; grid on;
% Lines to Earth
plot([0 rMars2E_max(1)], [0 rMars2E_max(2)], 'm-', 'LineWidth', 1.5, 'DisplayName', 'Line to Earth @ MAX');
plot([0 rMars2E_min(1)], [0 rMars2E_min(2)], 'c-', 'LineWidth', 1.5, 'DisplayName', 'Line to Earth @ MIN');
% Earth markers
plot(rMars2E_max(1), rMars2E_max(2), 'd', 'MarkerSize', 10, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'DisplayName', 'Earth @ MAX');
plot(rMars2E_min(1), rMars2E_min(2), 'd', 'MarkerSize', 10, 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c', 'DisplayName', 'Earth @ MIN');
legend('Location', 'best');
title('Spacecraft-to-Earth Vectors');
xlabel('X [km]'); ylabel('Y [km]');
axis equal;
% annotate distances
text(rMars2E_max(1)*0.6, rMars2E_max(2)*0.6, sprintf('%.1e km', maxInfo.range), 'Color', 'm', 'FontWeight', 'bold');
text(rMars2E_min(1)*0.6, rMars2E_min(2)*0.6, sprintf('%.1e km', minInfo.range), 'Color', 'c', 'FontWeight', 'bold');

sgtitle('Mars Orbiter & Earth Range Extremes');