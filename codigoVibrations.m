%% 
%% 
S = load('sensorlog_20260928_152835.mat');
A = S.Acceleration; G = S.AngularVelocity;
tA = seconds(A.Timestamp - A.Timestamp(1));
tG = seconds(G.Timestamp - G.Timestamp(1));
figure

subplot(2,1,1); plot(tA,[A.X A.Y A.Z]); grid on ;
legend('X','Y','Z'); ylabel('Aceleracao (m/s^2)');

subplot(2,1,2) ; plot(tG,[G.X G.Y G.Z]) ; grid on ;
legend('X','Y','Z'); ylabel('Vel. angular (rad/s)');
% angular speed of interest is given in the direction Z (orthogonal axis
% exiting the plan of the pendulum)

xlabel('Tempo (s)');

% --- Filter (Ignore beginning and ending) ---
% Window without human interference 
tempo_min = 1.0;  % letting go
tempo_max = 20.0; % stopping

% Cut in acceleration
masqueA = (tA >= tempo_min) & (tA <= tempo_max);

tA = tA(masqueA) - tempo_min; % Restarting vector time (zero)
A = A(masqueA, :); 

% Cut in angular speed
masqueG = (tG >= tempo_min) & (tG <= tempo_max);
tG = tG(masqueG) - tempo_min; % Restarting (useful time)
G = G(masqueG, :);        

% 1. Align the timelines so you can plot them together
% This interpolates the Z angular velocity to match the accelerometer timestamps
GZ_aligned = interp1(tG, G.Z, tA, 'linear', 'extrap');

% 2. True Radial Axis = Whichever axis (X or Y) sits right at 9.81 
% when the pendulum is perfectly still = RADIAL axis
A_radial = A.Y; 

% 4. Comparing !
% Find the first index where time surpasses 16 seconds
idx_start = find(tA > 16, 1, 'first');

% Find the last index where time is 20 seconds at most
idx_end = find(tA <= 18, 1, 'last');

% Extract the radial acceleration values between these two indices
A_radial_segment = A_radial(idx_start:idx_end);

% Calculate the average value to serve as actual g
g_mesure = mean(A_radial_segment);



% 1. Numerically integrate angular velocity to get the angle theta (in radians)
% Assuming tA is your time vector and GZ_aligned is in rad/s
theta = cumtrapz(tA, GZ_aligned); 

% Remove any linear drift from integration if necessary using detrend
theta = detrend(theta); 

% 2. True Dynamic Theoretical Radial Acceleration
R = 0.37; % Approximate distance to smartphone center of mass
% Instead of adding flat g_mesure, multiply it by the cosine of the shifting angle
Acp_theoretical = (GZ_aligned.^2) * R + abs(g_mesure) * cos(theta);



% Display the results in the command window
fprintf('Extracted g value: %.4f m/s²\n', g_mesure);
figure
plot(tA, abs(A_radial), 'b', 'DisplayName', 'Aceleração radial medida'); 
hold on;
plot(tA, Acp_theoretical, 'r--', 'DisplayName', 'Teorica (\omega^2*R + g*cos(\theta))'); 
%the value of g must be added since the smartphone sensor has a non zero
%value reference (the gravitation acceleration at neutral position)
hold on;
yline(-g_mesure, '--g', 'LineWidth', 1,'DisplayName','Valor absoluto estimado da gravidade'); 
grid on;
legend;
ylabel('Aceleracao (m/s^2)');
xlabel('Tempo (s)');
title('Validação da Aceleração Centripeta');

% the correct acceleration component is given in the Z direction (negative
% sense), because of the value being close to 9.81 m/s². It fluctuates
% around 9.81 m/s² because the neutral position has already this value (the
% gravitational field acts during all the data acquisition). The
% fluctuation must happen because of the rotational movement (centripetal
% acceleration takes place) :
% acp = v²/R = w².R, with v the tangential speed (maximum at the lowest point of
% the trajectory) and R the radios (equivalent to the distance from the
% center of mass of the smartphone used to the rotational axis fixed -
% pivoting distance). We cannot plot (tG, (G.Z.^2)/(.37^2)) with the
% Acceleration using the same graph to see if the equation hits the same values of Ar using the
% angular speed on Z axis as entry values. This is because tG is not the
% same as tA (software was not able to sample both using only one time
% vector)

% At the lowest point Fcp = P -> m.v²/R = m.g -> w².R = g -> 
% 2.pi.f = sqrt(g/R) -> f = sqrt(g/R)/2.pi (the formula for the
% natural frequency of a simple pendulum). Hence, at the lowest point, we
% should have w = sqrt(g/R) = sqrt (9.81 m/s ² / 29.5cm) = 5.7666 . 1/s =
% 5.766 Hz. However, our system is not a simple pendulum, but a 
% concentrated mass with a bar as the connection to the pivoting center).

barLength = 38/100;
barBase = 3/100;

%A30 SM-A305GT series RQ8M604MZFB IMEI 354881103280520 smartphone info
%https://www.maiscelular.com.br/fichas-tecnicas/samsung/galaxy-a30/sm-a305gtds/
phoneHeight = 15.85/100;
phoneBase = 7.47/100;

b1 = phoneBase ; h1 = phoneHeight ; 
b2 = barBase ; h2 = barLength ; 

y1 = -((h2 - 1/100) - h1/2); y2 = - ((h2 - 1/100)/2); 
% the pivoting center is admitted as the zero coordinate, explaining the
% negative values.

% If I was the inertia moment of AREA
%A1 = b1*h1 ; A2 = b2*h2; 
%h = (A1*y1 + A2*y2)/(A1 + A2); % system's centroid 
%d1 = abs(y1 - h) ; d2 = abs(y2 - h) ; 
%I1 = (b1*h1^3)/12 ; I2 = (b2*h2^3)/12;
%It = (I1 + A1*d1^2) + (I2 + A2*d2^2); 

%For I as the inertia moment of MASS, because period's formula has mass in
%the divisor
mS = 0.165 ; mb = 0.100; g = abs(g_mesure) ;
% Perpendicular Axis Theorem (the pivot is located on the Z-axis, and not
% on the X or Y axes)
I1 = (mS* (b1^2 + h1^2))/12 ; I2 = (mb*(b2^2 + h2^2))/12;

% Parallel Axis Theorem to the pivot (0,0)
d1 = abs(y1) ; d2 = abs(y2);
It = (I1 + mS*d1^2) + (I2 + mb*d2^2); 

h = (y1*mS + y2*mb) / (mb+mS) ;
h = abs(h) ; 
To = 2*pi*sqrt(It/((mS+mb)*g*h)); % kg*m² / (kg*m/s²*m)
fprintf('Calculated frequency: %.3f Hz\n', 1/To);

% Smartphone/Software sampling frequencies in Hz
fsA = 1/median(diff(tA)); fsG = 1/median(diff(tG));
% For the acceleration (A) is 125 Hz and for the angular velocity (G), we
% have 47.6190 Hz. Software sampling was fixed at 50 Hz.


% Estimate dominant pendulum frequency from the Z-axis angular velocity
angularVelocityZ = detrend(G.Z);
[powerSpectrum, frequency] = periodogram(angularVelocityZ, [], [], fsG,'power');
%instead of the traditional FFT, the periodogram is needed for sampling
%periods that fluctuate over time
[~, peakIndex] = max(powerSpectrum(2:end));
dominantFrequency = frequency(peakIndex + 1);
fprintf('Dominant frequency for: %.3f Hz\n', dominantFrequency);

figure
plot(frequency, powerSpectrum); grid on;
ylabel('Power (rad/s)² ')
% angular speed of interest is given in the direction Z (orthogonal axis
% exiting the plan of the pendulum)

xlabel('Frequência (Hz)')

%% ========================================================================
%% PARTE COMPLEMENTAR: ANÁLISE DE CICLOS E VALIDAÇÃO EXPERIMENTAL
%% ========================================================================

% 1. Detecção automática de picos para velocidade angular (G.Z)
% Usamos uma distância mínima entre picos baseada na frequência dominante aproximada
minPeakDistance_G = 0.8 * (1 / dominantFrequency); 
%0.8 para evitar de contar picos falsos
[pksG, locsG] = findpeaks(G.Z, tG, 'MinPeakDistance', minPeakDistance_G);
%pksG em rad/s e locsG em s

% 2. Execução das contas para a Questão 3: Método dos k-ciclos (Giroscópio)

k = 5;

% Ensaio 1: Ciclos de 1 a 6
t_p0_1 = locsG(1);
t_pk_1 = locsG(1 + k); % Índice 6
T_med1 = (t_pk_1 - t_p0_1) / k;
f_exp1 = 1 / T_med1;

% Ensaio 2: Ciclos de 6 a 11 (Totalmente independente do grupo 1)
t_p0_2 = locsG(1 + k);   % Índice 6
t_pk_2 = locsG(1 + 2*k); % Índice 11
T_med2 = (t_pk_2 - t_p0_2) / k;
f_exp2 = 1 / T_med2;


% Cálculo da dispersão absoluta entre ensaios manuais
dispersao_f = abs(f_exp1 - f_exp2);

% Print dos resultados da velocidade angular
fprintf('\n--- ANÁLISE DO GIROSCÓPIO (G_Z) ---\n');
fprintf('Ensaio 1 (Man): t0 = %.3fs, tk = %.3fs, T_med = %.3fs, f = %.3f Hz\n', t_p0_1, t_pk_1, T_med1, f_exp1);
fprintf('Ensaio 2 (Man): t0 = %.3fs, tk = %.3fs, T_med = %.3fs, f = %.3f Hz\n', t_p0_2, t_pk_2, T_med2, f_exp2);
fprintf('Dispersão absoluta experimental (f_exp1 - f_exp2): %.6f Hz\n', dispersao_f);


% 3. Estimativa independente através do sinal do acelerômetro (Questão 3)
% O acelerômetro radial deve apresentar o dobro de picos por ciclo
minPeakDistance_A = 0.4 * (1 / dominantFrequency);
[pksA, locsA] = findpeaks(abs(A_radial), tA, 'MinPeakDistance', minPeakDistance_A, 'MinPeakHeight', 9.5);

% Ensaio com o acelerômetro para k=10 picos (equivalente a 5 ciclos completos do pêndulo)
k_A = 10;
t_p0_A = locsA(1);
t_pk_A = locsA(1 + k_A);
T_med_A = (t_pk_A - t_p0_A) / (k_A / 2); % Dividido por (k_A/2) devido ao dobro de picos
f_exp_A = 1 / T_med_A;

fprintf('\n--- ANÁLISE INDEPENDENTE DO ACELERÔMETRO (A_Y) ---\n');
fprintf('Acelerômetro: t0 = %.3fs, tk = %.3fs, f_estimada = %.3f Hz (2x picos por ciclo)\n', t_p0_A, t_pk_A, f_exp_A);


% 4. Análise de Grande vs. Pequena Amplitude (Questão 4)
% Grande amplitude: medida nos primeiros ciclos capturados
T_grande = locsG(2) - locsG(1); 

% Pequena amplitude: medida nos últimos ciclos estáveis do registro
T_pequena = locsG(end) - locsG(end-1); 

% Diferença percentual em relação à pequena amplitude
dif_percentual_T = (abs(T_grande - T_pequena) / T_pequena) * 100;

fprintf('\n--- ANÁLISE DE AMPLITUDE (Questão 4) ---\n');
fprintf('Período em Grande Amplitude (Início): %.4f s\n', T_grande);
fprintf('Período em Pequena Amplitude (Fim): %.4f s\n', T_pequena);
fprintf('Diferença Percentual: %.3f %%\n', dif_percentual_T);


%% ========================================================================
%% VISUALIZAÇÃO GRÁFICA DOS PICOS ESCOLHIDOS
%% ========================================================================
figure
subplot(2,1,1)
plot(tG, G.Z, 'k'); hold on;
plot(locsG(1:7), pksG(1:7), 'ro', 'MarkerFaceColor', 'r');
grid on; title('Picos Selecionados no Giroscópio (G_Z)');
ylabel('Velocidade Angular (rad/s)'); xlabel('Tempo (s)');
legend('Sinal G_Z', 'Picos Marcados (k ciclos)');

subplot(2,1,2)
plot(tA, abs(A_radial), 'b'); hold on;
plot(locsA(1:12), pksA(1:12), 'mo', 'MarkerFaceColor', 'm');
grid on; title('Picos de Aceleração Radial (A_Y) - Frequência Dupla');
ylabel('Aceleração Radial (m/s²)'); xlabel('Tempo (s)');
legend('Sinal A_Y', 'Picos Inerciais (2 por ciclo)');



% --- DISPLAY QUANTITATIVE DATA SUMMARY ---
fprintf('\n================== RESULTS SUMMARY ==================\n');
fprintf('Dominant Frequency (Periodogram): %.3f Hz\n', dominantFrequency);
fprintf('Manual Gyro Group 1 Frequency:     %.3f Hz\n', f_exp1);
fprintf('Manual Gyro Group 2 Frequency:     %.3f Hz\n', f_exp2);
fprintf('Experimental Dispersion Score:    %.6f Hz\n', dispersao_f);
fprintf('Independent Accel Frequency:      %.3f Hz\n', f_exp_A);
fprintf('Large Amplitude Period (Start):   %.4f s\n', T_grande);
fprintf('Small Amplitude Period (End):     %.4f s\n', T_pequena);
fprintf('Non-linear Delta Variance:        %.3f %%\n', dif_percentual_T);