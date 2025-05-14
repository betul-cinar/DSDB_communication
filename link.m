clc; clear;
close all;

% Constants
k = 1.38e-23;              % Boltzmann's constant (J/K)
c = 3e8;                   % Speed of light (m/s)

% Mission Parameters
Pt_dBW = 27.78;               % Transmit power in dBW (20 W)

Gt_dBi = 56.6;               % Transmit antenna gain in dBi
Gr_dBi = 79.34;               % Receive antenna gain in dBi (DSN-36)
Lt_dB = 1.5;               % Transmit line loss in dB
Lp_dB = 1.0;               % Pointing loss in dB
La_dB = 2.0;               % Atmospheric loss in dB
Ts = 424;                  % System noise temperature in Kelvin
%R = 5.81e6;            % Data rate in bits per second
R = 2e6; 
f = 32e9;                  % Frequency (Hz)
lambda = c / f;           % Wavelength (m)
d = 394424035.6e3;                 % Distance (m) Mars-Earth
nu = 55/100;


% Free space path loss
Lfs_dB = 20*log10(4*pi*d*f/c);

% Total losses
L_total_dB = Lt_dB + Lp_dB + La_dB + Lfs_dB;

% EIRP
EIRP_dB = Pt_dBW + Gt_dBi - Lt_dB;

% Received power (in dBW)
Pr_dBW = EIRP_dB + Gr_dBi - L_total_dB;

% Convert received power to Watts
Pr_W = 10^(Pr_dBW/10);

% Compute Noise Power
No_W = k * Ts;
No_dBW_Hz = 10*log10(No_W);
N_dBW = No_dBW_Hz + 10*log10(R);  % Total noise over bandwidth

% Eb/No
EbNo_dB = Pr_dBW - N_dBW;

% Output
fprintf('--- Link Budget Results ---\n');
fprintf('EIRP             : %.2f dBW\n', EIRP_dB);
fprintf('Free Space Loss  : %.2f dB\n', Lfs_dB);
fprintf('Total Loss       : %.2f dB\n', L_total_dB);
fprintf('Received Power   : %.2f dBW\n', Pr_dBW);
fprintf('Noise Power      : %.2f dBW\n', N_dBW);
fprintf('Eb/No            : %.2f dB\n', EbNo_dB);
