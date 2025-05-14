clc; clear;

%% Constants
c = 3e8;                        % Speed of light [m/s]
h = 6.626e-34;                  % Planck's constant [J.s]
lambda = 1550e-9;               % Wavelength [m]

%% DSOC System Parameters
Pt_W = 4;                       % Transmit power [W]
D_tx = 0.22;                    % Transmitter aperture diameter [m]
D_rx = 5.1;                     % Receiver aperture diameter [m]
eta_tx = 0.8;                   % Transmitter optical efficiency
eta_rx = 0.8;                   % Receiver optical efficiency

% Distance
d_km = 394424035.6;             % Mars-Earth worst-case [km]
d_m = d_km * 1e3;               % Convert to meters

%% Optical Link Equation
% Pr = Pt * (Dt * Dr / (4 * lambda * d))^2 * eta_tx * eta_rx
gain_term = (D_tx * D_rx / (4 * lambda * d_m))^2;
Pr_W = Pt_W * gain_term * eta_tx * eta_rx;

%% Photon Flux
photons_per_sec = Pr_W * lambda / (h * c);

% Data rate
Rb = 10e6;  % 10 Mbps
photons_per_bit = photons_per_sec / Rb;

%% Display Results
fprintf('\n--- Optical Link Budget (Direct Model, 2.636 AU, 10 Mbps) ---\n');
fprintf('Distance to Earth         : %.2f AU (%.2e m)\n', d_km / 1.496e8, d_m);
fprintf('Transmit Power            : %.2f W\n', Pt_W);
fprintf('Gain Term                 : %.2e\n', gain_term);
fprintf('Received Power            : %.2e W\n', Pr_W);
fprintf('Photon Flux               : %.2e photons/sec\n', photons_per_sec);
fprintf('Data Rate                 : %.2f Mbps\n', Rb / 1e6);
fprintf('Photons per Bit           : %.2f\n', photons_per_bit);
fprintf('------------------------------------------------------------\n');
