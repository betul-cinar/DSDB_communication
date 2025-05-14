%===============================================================================
% Mars Orbiter -> DSN Line-of-Sight Analysis with Goldstone, Madrid, Canberra
%===============================================================================

clear; clc; close all

%% 1) Orbital & physical parameters
R_mars   = 3396.2;       % Mars radius [km]
alt      = 300;          % Orbiter altitude above Mars surface [km]
r_orb    = R_mars + alt; % Orbital radius [km]
mu_mars  = 4.282837e4;   % Mars gravitational parameter [km^3/s^2]

% Circular, equatorial orbit
incl     = 92;    % inclination [deg]
RAAN     = 0;     % RA of ascending node [deg]
argp     = 0;     % argument of periapsis [deg]
f0       = 0;     % initial true anomaly [deg]

%% 2) DSN ground station definitions (geodetic)
% Coordinates: [latitude (deg), longitude (deg, east+), altitude (km)]
stations = {
    'Canberra',  -35.4016,   148.9819, 0.7;   % DSS-43, Australia
    'Goldstone',  35.2472,  -116.7933, 1.0;   % DSS-14, CA, USA
    'Madrid',     40.4314,    -4.2486, 0.6;   % DSS-63, Spain
};
numSta = size(stations,1);

% WGS84 Earth ellipsoid constants
a_e = 6378.137;             % Equatorial radius [km]
f_e = 1/298.257223563;      % Flattening
e2  = 2*f_e - f_e^2;

% Precompute station ECEF (km)
stationECEF = zeros(3,numSta);
for i = 1:numSta
    lat = deg2rad(stations{i,2});
    lon = deg2rad(stations{i,3});
    h   = stations{i,4};
    N   = a_e / sqrt(1 - e2*sin(lat)^2);
    x   = (N + h)*cos(lat)*cos(lon);
    y   = (N + h)*cos(lat)*sin(lon);
    z   = (N*(1-e2) + h)*sin(lat);
    stationECEF(:,i) = [x; y; z];
end

%% 3) Time vector (UTC)
startTime = datetime(2032,07,25,0,0,0);
endTime   = datetime(2032,07,29,23,59,59);
dt        = 600;            % time-step [s]
tUTC      = (startTime:seconds(dt):endTime)';
N         = numel(tUTC);

%% 4) Preallocate
losSta     = false(numSta, N);  % line-of-sight matrix
theta      = zeros(N,1);        % orbiter true anomaly [rad]
rOrbitECI  = zeros(3,N);

%% 5) Mars orbiter ephem & station LOS loop
T_orb = 2*pi * sqrt(r_orb^3 / mu_mars);
n     = 2*pi / T_orb;

for k = 1:N
    JD          = juliandate(tUTC(k));
    rSun2Earth  = planetEphemeris(JD,'Sun','Earth')';
    rSun2Mars   = planetEphemeris(JD,'Sun','Mars')';
    rMars2Earth = rSun2Earth - rSun2Mars;

    dt_sec      = seconds(tUTC(k) - tUTC(1));
    theta(k)    = deg2rad(f0) + n * dt_sec;
    r_pf        = [r_orb*cos(theta(k)); r_orb*sin(theta(k)); 0];
    Rz_RAAN     = [cosd(RAAN) -sind(RAAN) 0; sind(RAAN) cosd(RAAN) 0; 0 0 1];
    Rx_INC      = [1 0 0; 0 cosd(incl) -sind(incl); 0 sind(incl) cosd(incl)];
    Rz_ARGP     = [cosd(argp) -sind(argp) 0; sind(argp) cosd(argp) 0; 0 0 1];
    rOrbitECI(:,k) = Rz_RAAN * Rx_INC * Rz_ARGP * r_pf;

    rSun2SC = rSun2Mars + rOrbitECI(:,k);

    T_ut1 = (JD - 2451545.0)/36525;
    GMST = 280.46061837 + 360.98564736629*(JD - 2451545.0) ...
           + 0.000387933*T_ut1^2 - T_ut1^3/38710000;
    GMST = mod(GMST,360);
    gst  = deg2rad(GMST);
    R_e  = [ cos(gst) sin(gst) 0; -sin(gst) cos(gst) 0; 0 0 1];

    for i = 1:numSta
        rSta_ECI = R_e * stationECEF(:,i);
        v = rSun2SC - (rSun2Earth + rSta_ECI);
        losSta(i,k) = (v' * rSta_ECI) > 0;
    end
end

%% 6) Aggregate daily contact time [hours]
days       = dateshift(tUTC,'start','day');
uniqueDays = unique(days);
ndays      = numel(uniqueDays);
contact_h  = zeros(numSta, ndays);
for j = 1:ndays
    idx = days == uniqueDays(j);
    for i = 1:numSta
        contact_h(i,j) = sum(losSta(i,idx))*dt/3600;
    end
end

%% 7) Display results
T = array2table(contact_h, 'RowNames', stations(:,1), ...
    'VariableNames', cellstr(datestr(uniqueDays,'yyyy_mm_dd')));
disp('Daily DSN contact times [hours]');
disp(T);

%% 8) Plot daily bars
figure('Name','Daily DSN Contact Time');

dataStack = contact_h([3,2,1], :)';  % ndays×3
b = bar(uniqueDays, dataStack, 'stacked');

% bottom segment → mustard
b(3).FaceColor = [0.9290 0.6940 0.1250];  % mustard
% middle segment → blue
b(2).FaceColor = [0 0.4470 0.7410];       % blue
% top segment → orange
b(1).FaceColor = [0.8500 0.3250 0.0980];   % orange

legend({'Canberra','Madrid','Goldstone'}, 'Location','eastoutside');
xlabel('Date'); ylabel('Contact Time [h]');
title('Daily DSN Contact Time for Mars Orbiter');
grid on;
datetick('x','dd-mmm-yyyy','keepticks');

