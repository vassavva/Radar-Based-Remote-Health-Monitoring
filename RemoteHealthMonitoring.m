
clear; clc; close all;

%% Add Directory of Files
Directory = 'C:\Users\vassa\OneDrive\Έγγραφα\Dissertation\Human Trials\3\1\';


%% Set Radar Configuration using Log File
logFile = fullfile(Directory, 'logFile.txt');
fileID = fopen(logFile, 'r');
if fileID == -1
    error('Log file not found.');
end

% Creates array with each line of the file
lines = {};
while ~feof(fileID)
    lines{end+1} = fgetl(fileID);
end
fclose(fileID);


%Finds Slope, number of samples and sample Rate and Chirp repetition
%Interval
profileLine = lines{contains(lines, 'ProfileConfig')};
profileVals = sscanf(profileLine(strfind(profileLine,'ProfileConfig,')+length('ProfileConfig,'):end), '%f,');
freqSlopeConst = profileVals(8); 
slope = 76.91;    % manually adding because convertion was wrong 
K = slope * 1e12;

idleTime = profileVals(3) * 1e-8;  % convert to seconds
adcStartTime = profileVals(5) * 1e-8;  % convert to seconds
rampEndTime = profileVals(6) * 1e-8;  % convert to seconds
Tr = 51.99e-6;

numSamples = profileVals(10);
sampleRate = profileVals(11) * 1e3;  % convert ksps to Hz
fprintf('Sample Rate: %.2f Hz\n', sampleRate);
fprintf('Slope: %.4f Hz\n', K);
fprintf('numSamples: %d\n', numSamples);

%Finds number of Chirps and number of Frames and frame Periodicity
frameLine = lines{contains(lines, 'FrameConfig,0,')};
frameVals = sscanf(frameLine(strfind(frameLine,'FrameConfig,')+length('FrameConfig,'):end), '%f,');
numChirps = frameVals(3);
numFrames = frameVals(4);
framePeriodicity_ms = frameVals(5) * 5e-9 * 1e3;  % convert to ms
fs = 1000 / framePeriodicity_ms;


fprintf('numChirps: %d\n', numChirps);
fprintf('numFrames: %d\n', numFrames);
fprintf('Frame rate: %.2f fps\n', fs);

%Finds number of Rx
chanLine = lines{contains(lines, 'ChannelConfig')};
chanVals = sscanf(chanLine(strfind(chanLine,'ChannelConfig,')+length('ChannelConfig,'):end), '%f,');
rxMask = chanVals(1);
txMask = chanVals(2);
numRX = sum(dec2bin(rxMask) == '1');  % count active RX bits
numTX = sum(dec2bin(txMask) == '1');  % count active TX bits
fprintf('numRX: %d\n', numRX);
fprintf('numTX: %d\n', numTX);

fprintf('framePeriodicity_ms: %f\n', framePeriodicity_ms);
%fprintf('All profileVals:\n');
%disp(profileVals);
%fprintf('profileLine raw: %s\n', profileLine);
%% Read Data from File
datafile= fullfile(Directory,'adc_data_Raw_0.bin');
fileID = fopen(datafile, 'rb');
if fileID == -1
    error('File not found. Check path.');
end

raw_data = fread(fileID,'int16');
fclose(fileID);
raw_data =(reshape(raw_data,8, []));  % one column per lane
fprintf('Raw data read successfully\n'); %sends command when finished reading raw data;

for collumn= 1:5
    fprintf('%d ', raw_data(:,collumn));
    fprintf('\n');
end



%% Seperating I and Q data
%I is real/horizontal part
%Q is imaginary/vertical part
% Extract I and Q for each channel
Rx0_I = raw_data(1, :);
Rx1_I = raw_data(2, :);
Rx2_I = raw_data(3, :);
Rx3_I = raw_data(4, :);

Rx0_Q = raw_data(5, :);
Rx1_Q = raw_data(6, :);
Rx2_Q = raw_data(7, :);
Rx3_Q = raw_data(8, :);

% Combine into complex samples per channel
Rx0 = complex(Rx0_I, Rx0_Q);
Rx1 = complex(Rx1_I, Rx1_Q);
Rx2 = complex(Rx2_I, Rx2_Q);
Rx3 = complex(Rx3_I, Rx3_Q);

complex_data = [Rx0; Rx1; Rx2; Rx3];%Only need Rx0


fprintf('complex_data:\n');
disp(complex_data(:, 1:5)');
%fprintf('%d ', Rx0_Q(1:5));

%reshape data

data = reshape(Rx0,numSamples, numChirps, []);%reshape data into a 3D matrix
numFrames = size(data,3);
fprintf('3D matrix:\n');
disp(data(1:5,1,1));

%% General Parameters

f_min = 77e9;
c = physconst('LightSpeed');
lambda_max = c/f_min;

B = K * Tr;             

%% Range FFT for slow time matrix (M x N)
rangeFFT = fft(data, numSamples, 1);
range_slow_time = squeeze(rangeFFT(:, 1, :))';  
 

%% Revoming of DC components using NNLS method
%plot IQ Constellation before and after 
range_slow_time_raw = range_slow_time;

for n = 1:numSamples
    samples = range_slow_time(:, n); %allocating the value of each range bin / collumn

    a = [real(samples), imag(samples)]; % create a vector with seperate values of the real and imaginary parts of the sapmles the samples
    A = [ones(numFrames,1) , -2*a]; %A5 equation
    b = -sum(abs(a).^2,2); %A5 equation

    y =  A\b;    % faster way of inv(A'*A) * A' * b - eq. A6 from paper

    ci = y(2);
    cq = y(3);

    range_slow_time(:,n) = samples - (ci + 1i*cq);

end

%% IQ Constellation Plot 
figure('Name', 'IQ Constellation Correction');

% Before DC compensation
scatter(real(range_slow_time_raw(:)), imag(range_slow_time_raw(:)), 5, [0.5 0.7 1], 'filled');
hold on;

% After DC compensation
scatter(real(range_slow_time(:)), imag(range_slow_time(:)), 5, [0.8 0.2 0.8], 'filled');

% Origin
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% Shifted origin (mean of raw signal)
shifted_origin = mean(range_slow_time_raw(:));
plot(real(shifted_origin), imag(shifted_origin), 'p', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');

axis equal; grid on;
xlabel('I'); ylabel('Q');
title('IQ Constellation Correction');
legend('Before DC compensation', 'After DC compensation', 'Origin', 'Shifted origin', ...
    'Location', 'best');

%% Angle Extraction and phase unwrappping
phase = angle(range_slow_time);
phase_unwrapped = unwrap(phase , [] ,1 );


%% Removal of DC value for each collumn
phase_no_dc = phase_unwrapped- mean(phase_unwrapped,1);

%% R0 Detection
range_mag = mean(abs(range_slow_time).^2, 1);  % average over frames then chirps
[~, range_bin]  = max(range_mag(2:numSamples)); 

fb = (range_bin)*(sampleRate)/numSamples;
R0 = (fb * c) / (2* K);                        % distance in metres
fprintf('R0: %.2f m\n', R0);

fprintf('range_bin: %d, fb: %.2f Hz, R0: %.2f m\n', range_bin, fb, R0);

rangeResolution = c / ( 2*K * Tr);
range_axis = (0:floor(numSamples)-1) * rangeResolution;

fprintf('Range resolution: %.2f m\n', rangeResolution);
fprintf('Range per bin: %.2f m\n', rangeResolution);
maxRange = rangeResolution * numSamples;
fprintf('Max range: %.2f m\n', maxRange);

%Range slow time

figure('Name', 'Range-Slow Time Map');

imagesc(range_axis, (0:numFrames-1)/fs, 20*log10(abs(range_slow_time)));
xlabel('Range (m)'); ylabel('Slow time (s)');
title('Range-Slow time matrix');
colorbar;

figure('Name','Range Profile');

plot(range_axis, 20*log10(range_mag(1:floor(numSamples))));
xlabel('Range (m)'); ylabel('Power (dB)');
title('Power Magnitude- Range ');
xline(R0, 'r--', sprintf('Detected R0=%.2fm', R0));
grid on;


figure('Name', 'Unwrapped Phase Map');
imagesc(range_axis, (0:numFrames-1)/fs, phase_unwrapped);
xlabel('Range (m)'); ylabel('Slow time (s)');
title('Unwrapped phase (Rad)');
colorbar;


%% Vibration FFT
vibrationFFT = fft(phase_no_dc, [] ,1); 

%%Range_Vibration Map

freqAxis = (0:numFrames-1) * (fs / numFrames);  % frequency axis in Hz
freqAxis_bpm = freqAxis * 60;
freq_mask = freqAxis_bpm(1:floor(numFrames/2)) >= 6 & freqAxis_bpm(1:floor(numFrames/2)) <= 120;

figure('Name', 'Range-Vibration Map');
vib_db = 20*log10(abs(vibrationFFT(1:floor(numFrames/2), :)) + eps);
imagesc(range_axis, freqAxis_bpm(freq_mask), vib_db(freq_mask, :));
clim([max(vib_db(freq_mask,:), [], 'all')-40, max(vib_db(freq_mask,:), [], 'all')]);
xlabel('Range (m)'); ylabel('Vibration frequency (per min)');
title('Range-Vibration map');
colorbar;

%% Range Bin 
range_vitals = freqAxis >= 0.1 & freqAxis <= 2.0;

rv_mag = abs(vibrationFFT(range_vitals, :)).^2;   % get magnitude of each complex number
averagePower = mean(rv_mag, 1);  % add this line

max_range_m = 5.0;  % only search up to 5m to avoid wall
max_range_bin = round(2*max_range_m / rangeResolution);
averagePower_limited = averagePower;
averagePower_limited(max_range_bin:end) = 0;
[~, range_bin_rv] = max(averagePower_limited);

fprintf('range_bin_rv: %d\n', range_bin_rv);

figure('Name', 'Average Power in Vital Signs Band');
plot((0:numSamples-1) * rangeResolution/2, averagePower);
xlabel('Range (m)'); ylabel('Average Power');
title('Average power in vital signs band');
xline(R0, 'r--', sprintf('R0=%.2fm', R0));
grid on;

%% Autocorrelation to extract HR and RR

fs_slow = fs;

% CREATE BANDPASS FILTERS
% Respiration filter: 0.1-0.5 Hz (paper uses 0.1-0.5 Hz)
bp_rr = designfilt('bandpassiir','FilterOrder', 4, 'HalfPowerFrequency1', 0.1, 'HalfPowerFrequency2', 0.5, 'SampleRate', fs_slow);


% Heart filter: 0.85-2.0 Hz
bp_hr = designfilt('bandpassiir','FilterOrder', 4, 'HalfPowerFrequency1', 0.85, 'HalfPowerFrequency2', 2.0, 'SampleRate', fs_slow);

% Extract phase at selected range bin
phase_signal = phase_no_dc(:, range_bin_rv);

% Apply bandpass filters
phase_rr = filtfilt(bp_rr, phase_signal);  % Respiratory component
phase_hr = filtfilt(bp_hr, phase_signal);  % Cardiac component

%% Autocorrelation for RR
[acf_rr, lags_rr] = xcorr(phase_rr, 'coeff');
acf_rr  = acf_rr(lags_rr >= 0);   % Keep only positive lags
lags_rr = lags_rr(lags_rr >= 0);

% PSD via Wiener-Khinchin theorem
PSD_rr       = abs(fft(acf_rr)).^2;
freq_axis_rr = (0:length(PSD_rr)-1) * (fs_slow / length(PSD_rr));

% Find RR peak 
rr_band = freq_axis_rr >= 0.1 & freq_axis_rr <= 0.5;
[~, idx_rr] = max(PSD_rr(:) .* rr_band(:));
RR_est = freq_axis_rr(idx_rr) * 60;


%% Autocorrelation for HR
[acf_hr, lags_hr] = xcorr(phase_hr, 'coeff');
acf_hr  = acf_hr(lags_hr >= 0);   % Keep only positive lags
lags_hr = lags_hr(lags_hr >= 0);

% PSD via Wiener-Khinchin theorem

PSD_hr       = abs(fft(acf_hr)).^2;
freq_axis_hr = (0:length(PSD_hr)-1) * (fs_slow / length(PSD_hr));

% Find HR peak 
hr_band = freq_axis_hr >= 0.833 & freq_axis_hr <= 2.0;
[~, idx_hr] = max(PSD_hr(:) .* hr_band(:));
HR_est = freq_axis_hr(idx_hr) * 60;

figure('Name', 'Vital Signs Waveforms');
time_axis = (0:numFrames-1) / fs;

subplot(2,1,1);
plot(time_axis, phase_rr);
xlabel('Time (s)'); ylabel('Displacement (rad)');
title('Respiration wave');
grid on;

subplot(2,1,2);
plot(time_axis, phase_hr);
xlabel('Time (s)'); ylabel('Displacement (rad)');
title('Heart wave');
grid on;



fprintf('RR: %.1f bpm, HR: %.1f bpm\n', RR_est, HR_est);

reference_RR = input('Enter reference RR (bpm): ');
reference_HR = input('Enter reference HR (bpm): ');

fprintf('Ref    RR: %.1f bpm\n', reference_RR);
fprintf('Ref    HR: %.1f bpm\n', reference_HR);

fprintf('Error  RR: %.1f%%\n', abs(RR_est - reference_RR) / reference_RR * 100);
fprintf('Error  HR: %.1f%%\n', abs(HR_est - reference_HR) / reference_HR * 100);

subplot(2,1,1);
plot(freq_axis_rr * 60, PSD_rr);
xlabel('Frequency (bpm)'); ylabel('Power');
title('RR PSD'); grid on;
xlim([6 60]);
xline(RR_est, 'r--', sprintf('RR=%.1f bpm', RR_est));
xline(reference_RR , 'g--', sprintf('Ref=%.1f bpm', reference_RR));

subplot(2,1,2);
plot(freq_axis_hr * 60, PSD_hr);
xlabel('Frequency (bpm)'); ylabel('Power');
title('HR PSD'); grid on;
xlim([48 120]);
xline(HR_est, 'r--', sprintf('HR=%.1f bpm', HR_est));
xline(reference_HR , 'g--', sprintf('Ref=%.1f bpm', reference_HR));

