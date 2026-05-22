%% ============================================================
% DAMPED FORCED OSCILLATOR
% Comparison between:
% 1) Convolution integral method
% 2) Direct numerical integration (ode45)
% 3) Frequency-domain method (FFT)
%
% Excitation: trapezoidal road bump of duration T
%
% Outputs:
% - time history of excitation
% - system responses
% - comparison between methods
% - frequency spectra
%% ============================================================

clear; close all; clc;

%% ------------------------------------------------------------
% 1. SYSTEM PARAMETERS (TOYOTA CELICA ST202 - 1994)
%% ------------------------------------------------------------

m = 400;              % mass [kg]
g = 9.81;             % gravity [m/s^2] (not used directly here)
k = 24516.626;        % stiffness [N/m]
c = 2985.65;          % viscous damping [N s/m]

zeta = (c/2) / sqrt(k*m);   % damping ratio [-]
wn = sqrt(k/m);              % natural angular frequency [rad/s]
wd = wn * sqrt(1 - zeta^2);  % damped natural frequency [rad/s]
fn = wn / (2*pi);           % natural frequency [Hz]

%% ------------------------------------------------------------
% 2. TIME DISCRETIZATION
%% ------------------------------------------------------------

dt = 0.001;        % time step [s]
Tmax = 10.0;       % total simulation time [s]

t = 0:dt:Tmax;     % time vector
N = length(t);

fs = 1/dt;         % sampling frequency [Hz]

%% ------------------------------------------------------------
% 3. EXCITATION: TRAPEZOIDAL ROAD BUMP
%% ------------------------------------------------------------

% Bump parameters
L = 0.9;        % bump length [m]
h = 0.03;       % bump height [m]
v_kmh = 10;     % vehicle speed [km/h]

% Unit conversion
v = v_kmh / 3.6;  % speed [m/s]

% Time characteristics
T = L / v;            % total excitation duration [s]
t_ramp = 0.25 * T;    % ramp-up / ramp-down time [s]
t_const = 0.75 * T;   % constant plateau duration reference

% Preallocation
y = zeros(size(t));
y_dot = zeros(size(t));

%% ------------------------------------------------------------
% Trapezoidal bump profile
%% ------------------------------------------------------------

for i = 1:length(t)
    t_i = t(i);

    if t_i >= 0 && t_i < t_ramp
        % ramp up
        y(i) = (h / t_ramp) * t_i;
        y_dot(i) = h / t_ramp;

    elseif t_i >= t_ramp && t_i < t_const
        % constant plateau
        y(i) = h;
        y_dot(i) = 0;

    elseif t_i >= t_const && t_i <= T
        % ramp down
        y(i) = (h / t_ramp) * (T - t_i);
        y_dot(i) = -h / t_ramp;

    else
        % outside bump
        y(i) = 0;
        y_dot(i) = 0;
    end
end

%% ------------------------------------------------------------
% Equivalent forcing: F(t) = c*y_dot + k*y
%% ------------------------------------------------------------

F = c .* y_dot + k .* y;

%% ------------------------------------------------------------
% Plot: bump profile
%% ------------------------------------------------------------

figure('Color','w','Name','Road Bump Profile');
plot(t, y, 'LineWidth', 2);
grid on;
xlim([0 2.5]);

xlabel('Time [s]');
ylabel('Vertical displacement [m]');
title('Trapezoidal road bump profile');

%% ------------------------------------------------------------
% 4. RESPONSE USING CONVOLUTION INTEGRAL
%% ------------------------------------------------------------

% Impulse response of SDOF system
h = (1/(m*wd)) * exp(-zeta*wn*t) .* sin(wd*t);

% Linear convolution (numerical approximation of integral)
x_conv = conv(F, h) * dt;
x_conv = x_conv(1:N);

v_conv = gradient(x_conv, dt);
a_conv = gradient(v_conv, dt);

%% ------------------------------------------------------------
% 5. RESPONSE USING ODE45
%% ------------------------------------------------------------

% Interpolated forcing function
Ffun = @(tt) interp1(t, F, tt, 'linear', 0);

% State-space ODE
odefun = @(tt, y) [
    y(2);
    (Ffun(tt) - c*y(2) - k*y(1)) / m
];

y0 = [0; 0];

[t_ode, y_ode] = ode45(odefun, [0 Tmax], y0);

% Resampling on uniform grid
x_ode = interp1(t_ode, y_ode(:,1), t, 'linear', 0);
v_ode = interp1(t_ode, y_ode(:,2), t, 'linear', 0);
a_ode = gradient(v_ode, dt);

%% ------------------------------------------------------------
% 6. FREQUENCY DOMAIN SOLUTION (FFT)
%% ------------------------------------------------------------

Nfft = 2^nextpow2(2*N - 1);

F_pad = zeros(1, Nfft);
F_pad(1:N) = F;

f_fft = (0:Nfft-1) * (fs / Nfft);
w_fft = 2*pi*f_fft;

% Frequency response function
H_fft = 1 ./ (k - m*w_fft.^2 + 1i*c*w_fft);

F_fft_full = fft(F_pad);
X_fft_full = H_fft .* F_fft_full;

x_fft_pad = real(ifft(X_fft_full)) * 2;
x_fft = x_fft_pad(1:N);

v_fft = gradient(x_fft, dt);
a_fft = gradient(v_fft, dt);

%% ------------------------------------------------------------
% 7. ERROR ANALYSIS
%% ------------------------------------------------------------

err_conv_ode = x_conv - x_ode;
err_fft_conv = x_fft - x_conv;
err_fft_ode  = x_fft - x_ode;

errmax_conv_ode = max(abs(err_conv_ode));
errmax_fft_conv = max(abs(err_fft_conv));
errmax_fft_ode  = max(abs(err_fft_ode));

fprintf('\n====================================================\n');
fprintf('                SYSTEM SUMMARY\n');
fprintf('====================================================\n');
fprintf('m      = %.2f kg\n', m);
fprintf('fn     = %.4f Hz\n', fn);
fprintf('wn     = %.4f rad/s\n', wn);
fprintf('zeta   = %.4f\n', zeta);
fprintf('k      = %.2f N/m\n', k);
fprintf('c      = %.2f N s/m\n', c);
fprintf('wd     = %.4f rad/s\n', wd);
fprintf('v      = %.2f m/s\n', v);

fprintf('\n====================================================\n');
fprintf('                 ERROR SUMMARY\n');
fprintf('====================================================\n');
fprintf('|x_conv - x_ode| max = %.6e m\n', errmax_conv_ode);
fprintf('|x_fft  - x_conv| max = %.6e m\n', errmax_fft_conv);
fprintf('|x_fft  - x_ode | max = %.6e m\n', errmax_fft_ode);


%% ------------------------------------------------------------
% 8. SPECTRAL ANALYSIS
%% ------------------------------------------------------------

f = (0:N-1)/(N*dt);
nHalf = floor(N/2) + 1;
f_plot = f(1:nHalf);

% Magnitude spectrum normalization helper
scaleSpectrum = @(X) (abs(X(1:nHalf))/N);

Fmag = scaleSpectrum(fft(F));
Fmag(2:end-1) = 2*Fmag(2:end-1);

Xc_mag = scaleSpectrum(fft(x_conv));
Xc_mag(2:end-1) = 2*Xc_mag(2:end-1);

Xf_mag = scaleSpectrum(fft(x_fft));
Xf_mag(2:end-1) = 2*Xf_mag(2:end-1);

Xo_mag = scaleSpectrum(fft(x_ode));
Xo_mag(2:end-1) = 2*Xo_mag(2:end-1);

%% ------------------------------------------------------------
% 9. PLOT: EXCITATION FORCE
%% ------------------------------------------------------------

figure('Color','w','Name','Forcing Function');
plot(t, F, 'LineWidth', 1.6);
grid on;
xlim([0 2.5]);

xlabel('Time [s]');
ylabel('Force [N]');
title('Excitation: Trapezoidal bump');

%% ------------------------------------------------------------
% 10. PLOT: CONVOLUTION RESPONSE
%% ------------------------------------------------------------

figure('Color','w','Name','Convolution Response');

plot(t, x_conv, 'LineWidth', 1.6);

grid on;
xlabel('Time [s]');
ylabel('x_{conv}(t) [m]');
title('Response using convolution integral method');

%% ------------------------------------------------------------
% 11. PLOT: FFT RESPONSE
%% ------------------------------------------------------------

figure('Color','w','Name','FFT Response');

plot(t, x_fft, 'LineWidth', 1.6);

grid on;
xlabel('Time [s]');
ylabel('x_{FFT}(t) [m]');
title('Response using frequency-domain (FFT) method');

%% ------------------------------------------------------------
% 12. PLOT: ODE45 RESPONSE
%% ------------------------------------------------------------

figure('Color','w','Name','ODE45 Response');

plot(t, x_ode, 'LineWidth', 1.6);

grid on;
xlabel('Time [s]');
ylabel('x_{ode45}(t) [m]');
title('Response using ODE45 numerical integration');

%% ------------------------------------------------------------
% 13. TIME HISTORY COMPARISON
%% ------------------------------------------------------------

figure('Color','w','Name','Time History Comparison');

plot(t, x_conv, 'LineWidth', 1.8); hold on;
plot(t, x_fft, '--', 'LineWidth', 1.5);
plot(t, x_ode, ':', 'LineWidth', 1.5);

grid on;

xlabel('Time [s]');
ylabel('x(t) [m]');
title('Comparison of system responses in time domain');

legend({'Convolution integral', 'FFT method', 'ODE45'}, ...
       'Location', 'best');

%% ------------------------------------------------------------
% 14. TIME HISTORY COMPARISON (ZOOM)
%% ------------------------------------------------------------

figure('Color','w','Name','Time History Comparison - Zoom');

plot(t, x_conv, 'LineWidth', 1.8); hold on;
plot(t, x_fft, '--', 'LineWidth', 1.5);
plot(t, x_ode, ':', 'LineWidth', 1.5);

grid on;

xlim([0 2.5]);

xlabel('Time [s]');
ylabel('x(t) [m]');
title('Comparison of responses (initial transient zoom)');

legend({'Convolution integral', 'FFT method', 'ODE45'}, ...
       'Location', 'best');

%% ------------------------------------------------------------
% 15. FORCE SPECTRUM
%% ------------------------------------------------------------

figure('Color','w','Name','Force Spectrum');

plot(f_plot, Fmag, 'LineWidth', 1.6);

grid on;

xlim([0 100]);

xlabel('Frequency [Hz]');
ylabel('|F(f)|');
title('Spectrum of excitation force');

%% ------------------------------------------------------------
% 16. RESPONSE SPECTRUM (CONVOLUTION METHOD)
%% ------------------------------------------------------------

figure('Color','w','Name','Response Spectrum - Convolution');

plot(f_plot, Xc_mag, 'LineWidth', 1.6);

grid on;

xlim([0 10]);

xlabel('Frequency [Hz]');
ylabel('|X_{conv}(f)|');
title('Response spectrum (convolution method)');

%% ------------------------------------------------------------
% 17. RESPONSE SPECTRUM (FFT METHOD)
%% ------------------------------------------------------------

figure('Color','w','Name','Response Spectrum - FFT');

plot(f_plot, Xf_mag, 'LineWidth', 1.6);

grid on;

xlim([0 10]);

xlabel('Frequency [Hz]');
ylabel('|X_{FFT}(f)|');
title('Response spectrum (FFT method)');

%% ------------------------------------------------------------
% 18. RESPONSE SPECTRUM (ODE45 METHOD)
%% ------------------------------------------------------------

figure('Color','w','Name','Response Spectrum - ODE45');

plot(f_plot, Xo_mag, 'LineWidth', 1.6);

grid on;

xlim([0 10]);

xlabel('Frequency [Hz]');
ylabel('|X_{ode45}(f)|');
title('Response spectrum (ODE45 method)');

%% ------------------------------------------------------------
% 19. SPECTRAL COMPARISON
%% ------------------------------------------------------------

figure('Color','w','Name','Spectral Comparison');

plot(f_plot, Xc_mag, 'LineWidth', 1.8); hold on;
plot(f_plot, Xf_mag, '--', 'LineWidth', 1.5);
plot(f_plot, Xo_mag, ':', 'LineWidth', 1.5);

grid on;

xlim([0 10]);

xlabel('Frequency [Hz]');
ylabel('|X(f)|');
title('Comparison of response spectra');

legend({'Convolution integral', 'FFT method', 'ODE45'}, ...
       'Location', 'best');

%% ------------------------------------------------------------
% 20. FFT vs CONVOLUTION ERROR (TIME DOMAIN)
%% ------------------------------------------------------------

figure('Color','w','Name','FFT vs Convolution Error');

plot(t, err_fft_conv, 'LineWidth', 1.6);

grid on;

xlim([0 2.5]);

xlabel('Time [s]');
ylabel('x_{FFT}(t) - x_{conv}(t) [m]');

title(sprintf('Difference between FFT and convolution methods (max error = %.3e m)', ...
      errmax_fft_conv));

%% ------------------------------------------------------------
% 21. CONVOLUTION vs ODE45 ERROR (TIME DOMAIN)
%% ------------------------------------------------------------

figure('Color','w','Name','Convolution vs ODE45 Error');

plot(t, err_conv_ode, 'LineWidth', 1.6);

grid on;

xlim([0 2.5]);

xlabel('Time [s]');
ylabel('x_{conv}(t) - x_{ode45}(t) [m]');

title(sprintf('Difference between convolution and ODE45 methods (max error = %.3e m)', ...
      errmax_conv_ode));

%% ------------------------------------------------------------
% 22. FFT vs ODE45 ERROR (TIME DOMAIN)
%% ------------------------------------------------------------

figure('Color','w','Name','FFT vs ODE45 Error');

plot(t, err_fft_ode, 'LineWidth', 1.6);

grid on;

xlim([0 2.5]);

xlabel('Time [s]');
ylabel('x_{FFT}(t) - x_{ode45}(t) [m]');

title(sprintf('Difference between FFT and ODE45 methods (max error = %.3e m)', ...
      errmax_fft_ode));

%% ------------------------------------------------------------
% 23. QUANTITATIVE METRICS (DISPLACEMENT RESPONSE)
%% ------------------------------------------------------------

% Preallocate results container (optional clarity)
results = struct();

%% --- MAXIMUM ABSOLUTE VALUE ---

[max_conv, idx_conv] = max(x_conv);
t_conv_max = t(idx_conv);

[max_fft, idx_fft] = max(x_fft);
t_fft_max = t(idx_fft);

[max_ode, idx_ode] = max(x_ode);
t_ode_max = t(idx_ode);

MaxAbs = max([max_conv, max_fft, max_ode]);

fprintf('\n====================================================\n');
fprintf('           MAXIMUM DISPLACEMENT\n');
fprintf('====================================================\n');
fprintf('Convolution : %.10f m at t = %.6f s\n', max_conv, t_conv_max);
fprintf('FFT         : %.10f m at t = %.6f s\n', max_fft, t_fft_max);
fprintf('ODE45       : %.10f m at t = %.6f s\n', max_ode, t_ode_max);
fprintf('----------------------------------------------------\n');
fprintf('Overall max : %.10f m\n', MaxAbs);

%% --- RMS VALUES ---

rms_conv = rms(x_conv);
rms_fft  = rms(x_fft);
rms_ode  = rms(x_ode);

fprintf('\n====================================================\n');
fprintf('                 RMS DISPLACEMENT\n');
fprintf('====================================================\n');
fprintf('Convolution : %.10f m\n', rms_conv);
fprintf('FFT         : %.10f m\n', rms_fft);
fprintf('ODE45       : %.10f m\n', rms_ode);

%% ------------------------------------------------------------
% 24. ACCELERATION RESPONSE (CONVOLUTION METHOD)
%% ------------------------------------------------------------

figure('Color','w','Name','Acceleration - Convolution');

plot(t, a_conv, 'LineWidth', 1.6);
grid on;
xlim([0 2.5]);

xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('Acceleration response - Convolution method');


%% ------------------------------------------------------------
% 25. ACCELERATION RESPONSE (FFT METHOD)
%% ------------------------------------------------------------

figure('Color','w','Name','Acceleration - FFT');

plot(t, a_fft, 'LineWidth', 1.6);
grid on;
xlim([0 2.5]);

xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('Acceleration response - FFT method');


%% ------------------------------------------------------------
% 26. ACCELERATION RESPONSE (ODE45 METHOD)
%% ------------------------------------------------------------

figure('Color','w','Name','Acceleration - ODE45');

plot(t, a_ode, 'LineWidth', 1.6);
grid on;
xlim([0 2.5]);

xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('Acceleration response - ODE45 method');


%% ------------------------------------------------------------
% 27. ACCELERATION TIME HISTORY COMPARISON
%% ------------------------------------------------------------

figure('Color','w','Name','Acceleration - Comparison');

plot(t, a_conv, 'LineWidth', 1.8); hold on;
plot(t, a_fft, '--', 'LineWidth', 1.5);
plot(t, a_ode, ':', 'LineWidth', 1.5);

grid on;
xlim([0 2.5]);

xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('Acceleration comparison');

legend({'Convolution', 'FFT', 'ODE45'}, 'Location', 'best');


%% ------------------------------------------------------------
% 28. ACCELERATION RMS VALUES
%% ------------------------------------------------------------

a_rms_conv = rms(a_conv);
a_rms_fft  = rms(a_fft);
a_rms_ode  = rms(a_ode);

fprintf('\n====================================================\n');
fprintf('              ACCELERATION RMS\n');
fprintf('====================================================\n');
fprintf('Convolution : %.10f m/s^2\n', a_rms_conv);
fprintf('FFT         : %.10f m/s^2\n', a_rms_fft);
fprintf('ODE45       : %.10f m/s^2\n', a_rms_ode);

%% ============================================================
% 29. CONTACT FORCE - CONVOLUTION METHOD
%% ============================================================

Ft_conv = c*(v_conv - y_dot) + k*(x_conv - y) + m*g;

figure('Color','w','Name','Contact Force - Convolution');

plot(t, Ft_conv, 'LineWidth', 1.6);
grid on;
xlim([0 5]);

xlabel('Time [s]');
ylabel('F_t^{conv}(t) [N]');
title('Contact force - Convolution method');

fprintf('\n====================================================\n');
fprintf('            CONTACT CONDITION CHECK\n');
fprintf('====================================================\n');

if all(Ft_conv >= 0)
    fprintf('Convolution : CONTACT OK\n');
else
    fprintf('Convolution : DETACHMENT DETECTED\n');
end


%% ============================================================
% 30. CONTACT FORCE - FFT METHOD
%% ============================================================

Ft_fft = c*(v_fft - y_dot) + k*(x_fft - y) + m*g;

figure('Color','w','Name','Contact Force - FFT');

plot(t, Ft_fft, 'LineWidth', 1.6);
grid on;
xlim([0 5]);

xlabel('Time [s]');
ylabel('F_t^{FFT}(t) [N]');
title('Contact force - FFT method');

if all(Ft_fft >= 0)
    fprintf('FFT         : CONTACT OK\n');
else
    fprintf('FFT         : DETACHMENT DETECTED\n');
end


%% ============================================================
% 31. CONTACT FORCE - ODE45 METHOD
%% ============================================================

Ft_ode = c*(v_ode - y_dot) + k*(x_ode - y) + m*g;

figure('Color','w','Name','Contact Force - ODE45');

plot(t, Ft_ode, 'LineWidth', 1.6);
grid on;
xlim([0 5]);

xlabel('Time [s]');
ylabel('F_t^{ode45}(t) [N]');
title('Contact force - ODE45 method');

if all(Ft_ode >= 0)
    fprintf('ODE45       : CONTACT OK\n');
else
    fprintf('ODE45       : DETACHMENT DETECTED\n');
end

fprintf('====================================================\n');