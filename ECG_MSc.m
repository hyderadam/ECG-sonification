%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ECG Analysis 
%%%%% Physionet -> MIT-BIH Arrhythmia Database -> 103m.mat
%%%%% 10/06/2018
%%%%% Adam Yassin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;  % Clear variables and functions from memory.
close all;  % Close figure.
clc;        % Clear command window.

%% Load ECG file %%
%%%%%%%%%%%%%%%%%%%

load('103m.mat'); % Load ECG data from MAT-file into workspace
raw_ECG = val(1,:);% raw_ECG parameter holds numerical values 
                    % of loaded ECG
% Data file (100m.mat) holds two signals. "val(1,:)" returns 
% the first of the two ECG signal from mat file (100m.mat).

% Convert raw ECG data to physical units by subtracting "base" 
% and dividing by "gain" ( base & gain are variables found in 
% info file of the ECG data)
base = 0; 
gain = 200;
% Formula for data conversion
ECG_signal = (raw_ECG - base) / gain;

Fs = 360;   % Sampling Hz -  found in info file of ecg data
Ts = 1/Fs;  % Sampling period
ECG_length = length(ECG_signal);        % Length of ECG data
t = 0:Ts:(ECG_length/Fs - Ts);          % Time axis (x-axis)
    
%% Cancellation of DC drift & Normalisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detrend signal
ECG_signal1 = detrend(ECG_signal);

% Normalise signal
ECG_signal1 = ECG_signal1 / max(abs(ECG_signal1)); 
 
%% Low-Pass Filter %%
%%%%%%%%%%%%%%%%%%%%%

% Transfer Function
%~ LPF ->  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2))

% Filter coefficients (Numerator)
b1 = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
% Filter coefficients (Denominator)
a1 = [1 -2 1];
% LPF Impulse Response
H_LPF = filter(b1,a1,[1 zeros(1,length(b1)-1)]);

% Returns vectors LPF and w1 containing the group delay, 
% and the frequencies
[LPF,w1] = grpdelay(b1,a1); 
                            
% Calculating LPF filter delay by averaging 
% the values from the group delay:
LPF_delay = round(nanmean(abs(LPF)));

% Applying LPF to the ECG signal
ECG_LPF = filter(b1,a1,ECG_signal1);
% Normalise signal
ECG_LPF = ECG_LPF / max(abs(ECG_LPF));

%% High-Pass Filter %%
%%%%%%%%%%%%%%%%%%%%%%

% Transfer Function
%~ HPF -> H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))

% Filter coefficients (Numerator)
b2 = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
% Filter coefficients (Denominator)
a2 = [1 1];
% HPF Impulse Response
H_HPF = filter(b2,a2,[1 zeros(1,length(b2)-1)]);

% Returns vectors HPF and w2 containing the group delay, 
% and the frequencies
[HPF,w2] = grpdelay(b2,a2);
                           
% Calculating HPF filter delay by averaging 
% the values from the group delay
HPF_delay = round(nanmean(abs(HPF)));

% Applying HPF to the ECG signal
ECG_HPF = filter(b2,a2,ECG_LPF);    
% Normalise signal
ECG_HPF = ECG_HPF / max(abs(ECG_HPF));

%% Derivative function %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Transfer Function
%~ derivative -> H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))

 % Filter coefficients (Numerator)
b3 = [1 2 0 -2 -1];
% Filter coefficients (Denominator)
a3 = 8*Ts;
% Derivative Function Impulse Response
H_derivative = filter(b3,a3,[1 zeros(1,length(b3)-1)]);

% Returns vectors derivative and w3 containing the group 
% delay and the frequencies
[derivative,w3] = grpdelay(b3,a3);

% Calculating derivative function delay by averaging 
% the values from the group delay
derivative_delay = round(nanmean(abs(derivative)));

% Applying derivative function to the ECG signal
ECG_derivative = filter(b3,a3,ECG_HPF);
% Normalise signal
ECG_derivative = ECG_derivative / max(abs(ECG_derivative));

%% Squaring function %%
%%%%%%%%%%%%%%%%%%%%%%%

% Multiplying every element in the modified ECG signal by 
% itself (Squaring the Signal)
ECG_SQR = ECG_derivative.^2; 
% Normalise signal
ECG_SQR = ECG_SQR / max(abs(ECG_SQR));

%% Moving window Function %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transfer Function
%~ window->Y(nt)=(1/N)[x(nT-(N - 1)T)+ x(nT - (N - 2)T)+...
%                                                     +x(nT)]

% Calculating moving window size -> ~ 150ms
windowSize = round(0.150*Fs);
% Filter coefficients (Numerator)
b4 = ones(1,windowSize)/windowSize;
% Filter coefficients (Denominator)
a4 = 1;
% Moving Window Function Impulse Response
H_window = filter(b4,a4,[1 zeros(1,length(b4)-1)]);

% Returns vectors window and w4 containing the group delay, 
% and the frequencies
[window,w4] = grpdelay(b4,a4);

% Calculating moving window function delay by averaging 
% the values from the group delay
window_delay = round(nanmean(abs(window)));

% Applying the average moving filter to the ECG signal
ECG_window = filter(b4,a4,ECG_SQR);
% Normalise signal
ECG_window = ECG_window / max(abs(ECG_window)); 

%% Peak Detection %%
%%%%%%%%%%%%%%%%%%%%

% Maximum value of the ECG window
max_window = max(ECG_window); 
% Average value of the ECG window
mean_window = mean(ECG_window );
% Calculating threshold value
threshold = max_window*mean_window; 

% Creates a vector (limit) of 0s and 1s, where 1s indicate 
% the values of the ECG signal that are greater than the threshold
limit =(ECG_window>(mean_window*max_window));

% Locates the position on the "limit" vector where 
% values go from 0 to 1.
QRS_start = find(diff([0 limit])==1);
% Locates the position on the "limit" vector where 
% values go from 1 to 0.
QRS_end = find(diff([limit 0])==-1);           
    
% Cancel Delay (because of pre-processing):
% ------------
% LPF delay 
% HPF delay
% Derivative filter delay
% window delay

% Readjusting values to aline with original ECG signal
QRS_start = QRS_start - (LPF_delay + HPF_delay + derivative_delay + window_delay);  
QRS_end = QRS_end - (LPF_delay + HPF_delay + derivative_delay + window_delay);     

% loop that checks if QRS_start vector has any negative 
% values. If it does, the loop will make that value = 0. 
% The reason for this is that the QRS_start vector may start 
% from a negative value, this is due to the realinement 
% (because of the offset cause by the filters). This occurs 
% because the first heart beat recorded by the ECG may not 
% start at the Q wave peak.
for idx = 1:numel(QRS_start)
    if QRS_start(idx) < 0
        QRS_start(idx) = 1;
    end
end

% Preparing peak variables
Q_val = zeros;
Q_loc = zeros;
R_val = zeros;
R_loc = zeros;
S_val = zeros;
S_loc = zeros;

% Peak Detection Loop

% - the loop is broken down into three parts:
% 1) Find the value and respective location of every R peak 
% in the ECG signal by obtaining the maximum value within 
% every limit box (from QRS_start to QRS_end). The R peak's 
% respective location is calculated by adding the discovered 
% location by its respective QRS_start location.
% 2) Find the value and respective location of every Q peak 
% in the ECG signal by obtaining the minimum value within 
% every limit box (from QRS_start to R_loc because the Q peak
% will always be before the R peak). The Q peak's respective 
% location is calculated by adding the discovered location by 
% its respective QRS_start location.
% 3) Find the value and respective location of every S peak 
% in the ECG signal by obtaining the minimum value within 
% every limit box (from R_loc to QRS_end because the S peak 
% will always be after the R peak). The S peak's respective 
% location is calculated by adding the discovered location by 
% its respective R peak location.
for i=1:length(QRS_start)
[R_val(i), R_loc(i)] = max(ECG_signal(QRS_start(i):QRS_end(i)));
R_loc(i) = R_loc(i) + QRS_start(i) - 1; % (-1 for loop offset)
 
[Q_val(i), Q_loc(i)] = min(ECG_signal(QRS_start(i):R_loc(i)));
Q_loc(i) = Q_loc(i) + QRS_start(i) - 1; % (-1 for loop offset)

[S_val(i), S_loc(i)] = min(ECG_signal(R_loc(i):QRS_end(i)));
S_loc(i) = S_loc(i) + R_loc(i) - 1;     % (-1 for loop offset)
end

%% Diagnosis %%
%%%%%%%%%%%%%%%

avg_QR_riseTime = mean(R_loc-Q_loc)/Fs;
avg_QR_riseLevel = mean(R_val-Q_val);
avg_RS_fallTime = mean(S_loc-R_loc)/Fs;
avg_RS_fallLevel = mean(R_val-S_val); 

disp('Average Q -> R rise time (seconds)');
disp(avg_QR_riseTime) % Average Rise time (in seconds)
disp('Average Q -> R rise level (mV)');
disp(avg_QR_riseLevel)   % Average Rise Level
disp('Average R -> S fall time (seconds)');
disp(avg_RS_fallTime) % Average Fall time (in seconds)
disp('Average R -> S fall level (mV)');
disp(avg_RS_fallLevel)  % Average Fall Level

% Determining Beats Per Minute (BPM)  
% Divide the beats (amount of R peaks) counted by the signal 
% duration (in minutes) to obtain BPM

% Find duration of ECG signal in seconds
duration_in_seconds = ECG_length/Fs;
% Find duration of ECG signal in minutes
duration_in_minutes = duration_in_seconds/60;
% Calculate BPM of signal.
BPM_avg = round(length(R_loc)/duration_in_minutes);

% Preparing variables for average QRS duration
qrs_avg = 0;
QRS_AVG = 0;

% Average QRS duration loop

% - calculates the distance between every Q peak and its 
% respective S peak (divided by Fs to obtain value in seconds) 
% in the ECG signal. The calculated values are added together 
% to then be divided by the number of R peaks, so as to obtain
% the average QRS duration of the signal.

for n = 1:length(S_loc)
    qrs_avg = ( S_loc(n) - Q_loc(n) ) / Fs;
    QRS_AVG = QRS_AVG  + qrs_avg;
end

% Multiplying by 1.69 to make up for lost duration 
% (half of Q and S wave are not accounted for)   
QRS_AVG = (QRS_AVG / length(R_loc))*1.69; 
                                                

% Loop that checks if average QRS duration is healthy
if QRS_AVG > 0.120
    disp('Averge QRS wave duration is greater than the healthy average (seconds)')
    disp(QRS_AVG)
elseif QRS_AVG < 0.08
    disp('Average QRS wave duration is smaller than the healthy average (seconds)')
    disp(QRS_AVG)
else
    disp('Average QRS wave duration is within a healthy boundary (seconds)')
    disp(QRS_AVG)
end

% Loop that checks if average Q peak value is healthy
if abs(mean(Q_val)) > abs(0.25 * mean(R_val))
    n1 = 1;
    disp('Averge Q peak amplitude is greater than the healthy average (mV)')
    disp(mean(Q_val))
elseif abs(mean(Q_val)) < abs(0.05 * mean(R_val))
    n1 = 1;
    disp('Averge Q peak amplitude is smaller than the healthy average (mV)')
    disp(mean(Q_val))
else
    disp('Average Q peak amplitude is within a healthy boundary (mV)')
    n1 = 0;
    disp(mean(Q_val))
end

% Loop that checks if average R peak value is healthy
if mean(R_val) > 1.7
    n2 = 1;
    disp('Averge R peak amplitude is greater than the healthy average (mV)')
    disp(mean(R_val))
elseif mean(R_val) < 0.18
    n2 = 1;
    disp('Averge R peak amplitude is smaller than the healthy average (mV)')
    disp(mean(R_val))
else
    disp('Average R peak amplitude is within a healthy boundary (mV)')
    n2 = 0;
    disp(mean(R_val))
end

% Loop that checks if average S peak value is healthy
if S_val < - 0.49
    n3 = 1;
    disp('Averge S peak amplitude is greater than the healthy average (mV)')
    disp(mean(S_val))
elseif mean(S_val) >= 0
    n3 = 1;
    disp('Averge S peak amplitude is smaller than the healthy average (mV)')
    disp(mean(S_val))
else
    disp('Average S peak amplitude is within a healthy boundary (mV)')
    n3 = 0;
    disp(mean(S_val))
end

% BPM - Intensity Alteration - Amplitude Modulation

% Takes the resultng sine wave from the frequency altering 
% step and modulates its amplitude sinusoidally. The rate of
% this modulation is dependent on the average beat per minute 
% of the ECG signal. If it is above the global average, the 
% signal is modulated at a higher frequency and the opposite 
% is true when the BPM average is lower than the 
% global average.
if BPM_avg >= 60 || BPM_avg <= 100
    x_f_morph = 2;
    disp('Average BPM is within a healthy boundary (BPM)')
elseif BPM_avg < 60
    x_f_morph = 1;
    disp('Average BPM indicates "Sinus Tachycardia" (BPM)')
elseif BPM_avg > 100
    x_f_morph = 4;
    disp('Average BPM indicates "Sinus bradycardia" (BPM)') 
end

% Dispay average BPM 
disp(BPM_avg);

% Display the amount of counted R peaks
disp('Number of detected R peaks')
disp(length(R_loc));

%% Data Sonification %%
%%%%%%%%%%%%%%%%%%%%%%%

% QRS - Frequency Alteration

% Produces a sine wave with a frequency that is based on
% the variable "morph" (which is a value that depends on how 
% close the average QRS duration of the signal is to a global 
% average). The close it is to the global average, the closer 
% the frequency will be to middle C.

% Standard QRS duration
QRS_AVG_True = 0.09;
% Fequency Altering paremeter
Hz = round((QRS_AVG - QRS_AVG_True)*2500);

x_fs = 44100;       % Sampling frequency
x_f1 = 262 + Hz;    % Sine wave frequency     
x_T = 1/x_fs;       % Period    
x_t = 0:x_T:3;      % Sound duration
x = sin(2*pi*x_f1*x_t);         % Sine wave
x = 0.9 .* (x / max(abs(x)));   % Normalise

% BPM - Intensity Alteration - Amplitude Modulation

morph = sin(2*pi*x_f_morph*x_t); % Sine wave of Amp modulation
morph =(morph / max(abs(morph)));% Normalise
x1 = 0.9 .* (x .* morph);        % Combination of both waves

% Q, R, and S values - Gibbs Phenomenon (ROUGHNESS)
nt = n1 + n2 + n3;  % Roughness parameter      

if nt == 0
        x2 = x1;
elseif nt == 1
        x2 = (x1 + ((1/3).*(sin(2*pi*x_f1*x_t*3)).*(morph)));
elseif nt == 2
        x2 = (x1 + ((1/3).*(sin(2*pi*x_f1*x_t*3)).*(morph)) + ((1/5).*(sin(2*pi*x_f1*x_t*5)).*(morph)));
elseif nt == 3
        x2 = (x1 + ((1/3).*(sin(2*pi*x_f1*x_t*3)).*(morph)) + ((1/5).*(sin(2*pi*x_f1*x_t*5)).*(morph)) + ((1/7).*(sin(2*pi*x_f1*x_t*7)).*(morph)));
end

% Normalise
x2 = 0.9 .* (x2 / max(abs(x2)));   

% Play Sound (diagnosis)
sound(x2,x_fs,16)

%% Plots %%
%%%%%%%%%%%

df = (Fs/2)/512;
xf = 0:df:Fs/2 - df;    %Creating x-axis.

% Maximize figure to screen size
set(figure(1), 'Position', get(0,'Screensize'));
figure(1)
axis tight;
subplot(4,1,1)
plot(t/60,ECG_signal)
title('ECG Signal');
xlabel('Time /min');
ylabel('Voltage /mV');

subplot(4,1,2)
plot(t(1:Fs*10),ECG_signal(1:Fs*10))
title('ECG signal')
xlabel('Time /s');
ylabel('Voltage /mV');
xlim([0 10])

subplot(2,1,2)
plot (t(1:Fs*5),ECG_signal(1:Fs*5),t(Q_loc) , Q_val, 'o',  t(R_loc) ,R_val , '^', t(S_loc) ,S_val, '*' );
xlim([0 5])
title('Analysed  ECG signal')
xlabel('Time /s');
ylabel('Voltage /mV');
legend('ECG','Q','R','S');

% Maximize figure to screen size
set(figure(2), 'Position', get(0,'Screensize'));
figure(2)
axis tight;
subplot(2,2,1)
plot(t(1:Fs*3),ECG_signal(1:Fs*3))
title('ECG Signal')
xlabel('Time /s');
ylabel('Voltage /mV');
xlim([0 3])

subplot(2,2,2)
plot(t(1:Fs*3),ECG_LPF(1:Fs*3))
title('LP filtered ECG Signal')
xlabel('Time /s');
ylabel('Voltage /mV');
xlim([0 3])

subplot(2,2,3)
zplane(b1,a1)
title('LPF poles and zeros')

subplot(2,2,4)
[mag_lpf,~] = freqz(b1,a1);
plot(xf,20*log10(abs(mag_lpf(1:512))));
title('LPF Magnitude Response')
xlabel('Frequency /Hz');
ylabel('Magnitude /dB');

% Maximize figure to screen size
set(figure(3), 'Position', get(0,'Screensize'));
figure(3)
axis tight;
subplot(2,2,1)
plot(t(1:Fs*3),ECG_LPF(1:Fs*3))
title('LP filtered ECG Signal')
xlabel('Time /s');
ylabel('Voltage /mV');
xlim([0 3])

subplot(2,2,2)
plot(t(1:Fs*3),ECG_HPF(1:Fs*3))
title('HP filtered ECG Signal')
xlabel('Time /s');
ylabel('Voltage /mV');
xlim([0 3])

subplot(2,2,3)
zplane(b2,a2)
title('HPF poles and zeros')

subplot(2,2,4)
[mag_lpf,~] = freqz(b2,a2);
plot(xf,20*log10(abs(mag_lpf(1:512))));
title('HPF Magnitude Response')
xlabel('Frequency /Hz');
ylabel('Magnitude /dB');

% Maximize figure to screen size
set(figure(4), 'Position', get(0,'Screensize'));
figure(4)
axis tight;
subplot(2,2,1)
plot(t(1:Fs*3),ECG_HPF(1:Fs*3))
title('HP filtered ECG Signal')
xlabel('Time /s');
ylabel('Voltage /mV');
xlim([0 3])

subplot(2,2,2)
plot(t(1:Fs*3),ECG_derivative(1:Fs*3))
title('Derivative filtered ECG Signal')
xlabel('Time /s');
ylabel('Voltage/s /mV/s');
xlim([0 3])

subplot(2,2,3)
zplane(b3,a3)
title('Derivative filter poles and zeros')

subplot(2,2,4)
[mag_lpf,~] = freqz(b3,a3);
plot(xf,20*log10(abs(mag_lpf(1:512))));
title('Derivative filter Magnitude Response')
xlabel('Frequency /Hz');
ylabel('Magnitude /dB');

% Maximize figure to screen size
set(figure(5), 'Position', get(0,'Screensize'));
figure(5)
axis tight;
subplot(2,2,1)
plot(t(1:Fs*3),ECG_SQR(1:Fs*3))
title('Squared derivative filtered ECG Signal')
xlabel('Time /s');
ylabel('Voltage/s /mV/s');
xlim([0 3])

subplot(2,2,2)
plot(t(1:Fs*3),ECG_window(1:Fs*3))
title('Moving average filter on ECG Signal')
xlabel('Time /s');
ylabel('Voltage/s /mV/s');
xlim([0 3])

subplot(2,2,3)
zplane(b4,a4)
title('Moving averge winodw poles and zeros')

subplot(2,2,4)
[mag_lpf,frq_lpf] = freqz(b4,a4);
plot(xf,20*log10(abs(mag_lpf(1:512))));
title('Average Moving Window Magnitude Response')
xlabel('Frequency /Hz');
ylabel('Magnitude /dB');

% Maximize figure to screen size
set(figure(6), 'Position', get(0,'Screensize'));
figure(6)
axis tight;

subplot(6,1,1)
plot(t(1:Fs*10),ECG_signal(1:Fs*10))
title('ECG signal')
xlabel('Time /s');
ylabel('Voltage /mV');

subplot(6,1,2)
plot(t(1:Fs*10),ECG_LPF(1:Fs*10))
title('Low-Pass')
xlabel('Time /s');
ylabel('Voltage /mV');

subplot(6,1,3)
plot(t(1:Fs*10),ECG_HPF(1:Fs*10))
title('High-Pass')
xlabel('Time /s');
ylabel('Voltage /mV');

subplot(6,1,4)
plot(t(1:Fs*10),ECG_derivative(1:Fs*10))
title('ECG signal')
xlabel('Time /s');
ylabel('Voltage/s /mV/s');

subplot(6,1,5)
plot(t(1:Fs*10),ECG_SQR(1:Fs*10))
title('ECG signal')
xlabel('Time /s');
ylabel('Voltage/s /mV/s');

subplot(6,1,6)
plot(t(1:Fs*10),ECG_window(1:Fs*10))
title('ECG signal')
xlabel('Time /s');
ylabel('Voltage/s /mV/s');
