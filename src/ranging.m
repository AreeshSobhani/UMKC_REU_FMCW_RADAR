clear;
clc;
close all;

%% parameters

c = 3e8;           % Speed of light
stratFreq = 77e9;  % Start frequency

Tr = 80e-6;          % Individual Chirp duration (also the period)
Samples = 512;       % Number of sampling points (configurable)
Fs = 1/Tr*Samples;   % Sampling rate

rangeBin = Samples;  % Number of range bins
chirps = 1024;        % Number of chirps (configurable)
dopplerBin = chirps; % Number of Doppler bins

slope = 29982e9;       % Chirp slope (configurable)
Bandwidth = slope * Tr;  % Effective bandwidth of transmitted signal
BandwidthValid = Samples/Fs*slope; % Actual signal bandwidth based on sampling
centerFreq = stratFreq + Bandwidth / 2; % Center frequency
lambda = c / centerFreq; % Wavelength

max_dist = (Fs*c)/(2*slope); % max usable radar distance

%channel number is product of rx and tx antennas -> "virtualAntenna"
txAntenna = ones(1,4); % Transmitting antenna configuration (configurable)
rxAntenna = ones(1,8); % Receiving antenna configuration (configurable)
txNum = length(txAntenna);
rxNum = length(rxAntenna);
virtualAntenna = length(txAntenna) * length(rxAntenna);

dz = lambda / 2; % Elevation spacing between receiving antennas (configurable)
dx = lambda / 2; % Horizontal spacing between receiving antennas (configurable)
% doaMethod = 2; % DOA estimation method: 1-DBF, 2-FFT, 3-Capon
target = [
    13 2 -30; % target1: range (m), speed (m/s), angle (degrees) , all these + num of targets is configurable
    5 -1 20
    20 4 10
    25 3 -15
    ];
targetNum = size(target,1);

%% generata signals
rawData = zeros(virtualAntenna,Samples,chirps); %pre-all 

t = 0:1/Fs:Tr-(1/Fs); % Time sequence for chirp sampling
for chirpId = 1:chirps %iterate for every chirp
    for txId = 1:txNum %iterate for each tx antenna
        St = exp((1i*2*pi)*(centerFreq*(t+((txNum-1)*chirps + chirpId)*Tr)+slope/2*t.^2)); % Transmitted signal
        for rxId = 1:rxNum %iterate for each rx antenna
            Sif = zeros(1,rangeBin);
            for targetId = 1:targetNum
                
                %3 TARGETS IN RADAR RANGE
                targetRange = target(targetId,1);
                targetSpeed = target(targetId,2);
                targetAngle = target(targetId,3);

                tau = 2 * (targetRange + targetSpeed * (txId - 1) * Tr) / c; % round trip delay
                fd = 2 * targetSpeed / lambda; % doppler frequency shift;
                wx = ((txId-1) * rxNum + rxId) / lambda * dx * sind(targetAngle); % phase shift?
                Sr = 10*exp((1i*2*pi)*((centerFreq-fd)*(t-tau+((txNum-1)*chirps + chirpId) * Tr)+slope/2*(t-tau).^2 + wx));  % Echo signal
                Sif = Sif + St .* conj(Sr); % IF signal 
            end
            rawData((txId-1) * rxNum + rxId,:,chirpId) = Sif; %store IF in raw data
        end
    end
end


%% Select One Chirp 
% plot the real part and imag part of that chirp
%separating real and imaginary components all signals

%extract single chirp IF signal
sing_chirp = rawData(1,:,1);

%extract magntiude and phase
chirp_mag = real(sing_chirp);
chirp_phase = imag(sing_chirp);

%plot chirp
figure(1);
title('Amplitude and Phase of a Single Chirp IF Signal');
xlabel('Time (s)');
ylabel('Amplitude (dB)');
hold on
plot(t,chirp_mag,'DisplayName','Magnitude');
plot(t,chirp_phase,'DisplayName','Phase');
hold off
legend

%% For All Chirps
% perform 1D FFT to see the range-Doppler Spectrum 
% using mesh() to plot 3D surface data with grid lines, and use view() to select top-down (2D-like) view. 

%pre-all
chirp_trans = zeros(Samples); 
chirp_trans_all = zeros(virtualAntenna,Samples,chirps);

%graph axes
speedRes = lambda / (2 * dopplerBin * Tr); %speed resolution
vel = (-dopplerBin/2:1:dopplerBin/2 - 1) * speedRes; %velocity axis
rangeRes = c / (2 * Bandwidth); % Range resolution
range = (0:rangeBin-1) * rangeRes; % Range axis

samp_data = squeeze(rawData(1,:,:)); %all samples for all chirps

wind1 = hanning(Samples); %hanning window function

%1d fft for every chirp
for chan_idx = 1:virtualAntenna
    for idx = 1:chirps
        chirp_trans(:,idx) = abs((fft(samp_data(:,idx).*wind1))); 
    end
    chirp_trans_all(chan_idx,:,:) = chirp_trans;
end

%range fft plot
figure(2);
imagesc(vel',range,db(squeeze(chirp_trans_all(1,:,:))));
title('Range-FFT');
xlabel('Velocity (m/s)');
ylabel('Range (m)');

%% 2D FFT Matrix
% get a 2D FFT Matrix for all chirps, this is a 3D matrix, like rawData(rangeId,angleIndex,dopplerId)

%pre-all
sig_sum = 0;
chirp_trans_2 = zeros(virtualAntenna,Samples,chirps);

%pre-all
chirps_trans_2 = zeros(virtualAntenna,Samples,chirps); 

%hanning window function in second dim
wind2 = hanning(chirps)';

%range-doppler calc and summing
for chan_num = 1:virtualAntenna
    samp_data_rang_dopp = squeeze(rawData(chan_num,:,:));
    chirp_trans_2(chan_num,:,:) = fftshift(fft2(samp_data_rang_dopp.*wind1.*wind2),2); %2d fft with window
    sig_sum = sig_sum + abs(chirp_trans_2(chan_num,:,:)); %sum magnitude of each peak
end

%% Range-Doppler 
% perform range-doppler to 2D FFT Matrix

%range-doppler plot
figure(3);
imagesc(vel, range, db(abs(squeeze(chirp_trans_2(1,:,:)))));
title('Range-Doppler Spectrum');
ylabel('Range (m)');
xlabel('Velocity (m/s)');

%% Range-Angle
% perform range-angle to 2D FFT Matrix

angleRes = lambda / (virtualAntenna * dx) * 180 / pi; % Angular resolution (degrees)
angleIndex = -(-virtualAntenna/2:1:virtualAntenna/2 - 1) * angleRes; %angle bins

wind3 = hanning(virtualAntenna);
wind4 = hanning(rangeBin)';

ang_sum = 0;
ran_ang_fft = zeros(virtualAntenna,Samples,chirps);

%range-angle fft calc
for chirp_num = 1:chirps
    samp_data_rang_ang = squeeze(rawData(:,:,chirp_num));
    anglefft = fftshift(fft(samp_data_rang_ang.*wind3.*wind4));
    ran_ang_fft(:,:,chirp_num) = fft(anglefft,[],2);
    ang_sum = ang_sum + abs(ran_ang_fft(:,:,chirp_num));
end

%range fft plot
figure(4);
imagesc(range, angleIndex, db(abs(ran_ang_fft(:,:,1))));
title('Range-Angle Spectrum');
xlabel('Range (m)');
ylabel('Angle of Attack (°)')

%% channel Accumulate
% sum signals from multiple channels (e.g., antennas, chirps, or pulses) to improve signal-to-noise ratio (SNR) 

%DONE WITHIN OG RANGE-DOPPLER LOOP

%% Generate CFAR Data 
% get a Range-Doppler Map

sig_sum = squeeze(sig_sum); %removing singleton dim

%summed range doppler plot
figure(5);
imagesc(vel, range, db(abs(sig_sum)));
title('Summed Range-Doppler Spectrum');
ylabel('Range (m)');
xlabel('Velocity (m/s)');

%% Peak Search 
% locate target peaks

%reference vals for cfar 
bias = 5; %mult factor for thres
guard = 1; %ignored cells next to COI
ref = 3; %cells used for thres
buffer = ref+guard; %buffer from COI
min_ind = buffer+1; %minimum usable val
max_ind = Samples-buffer; %maximum usable val

%pre-all
thres = zeros(1,max_ind);
range_detect = zeros(1,max_ind);
vel_detect = zeros(1,max_ind);
ang_detect = zeros(1,max_ind);

%cfar alg
for cfar_samp = min_ind:max_ind       
    %finding ref cells
    left_ref = mean(sig_sum(((cfar_samp-(ref+1)):(cfar_samp-(guard+1))),:),2); 
    right_ref = mean(sig_sum(((cfar_samp+(guard+1)):(cfar_samp+(ref+1))),:),2);
    ref_vec = cat(1,left_ref,right_ref);

    %pads unused indexes w 0
    if (length(ref_vec) == ref) | (length(ref_vec) == 2*ref)
        thres(cfar_samp) = mean(ref_vec)*bias;
    else
        thres(cfar_samp) = 0;
    end 
    
    %peak is a target if it passes the thres val at its index
    if mean(sig_sum(cfar_samp,:)) > thres(cfar_samp)
        range_detect(cfar_samp) = range(:,cfar_samp);
    end
end

%% DOA (Direction of Arrival)
% finally get the distance, angle, and speed of the targets 

range_vals = findpeaks(range_detect); %target range values

%pre-all
find_idx = zeros(1,length(range_vals));

%find peak vals for vel and AoA using range indexes
for det = 1:length(range_detect)
    if range_detect(det) ~= 0
        [~,max_y] = find(sig_sum == max(sig_sum(det,:)));
        vel_detect(det) = vel(max_y);
        [max_x,~] = find(ang_sum == max(ang_sum(:,det)));
        ang_detect(det) = angleIndex(max_x);
    end
end

%find indexes for vel and ang values to find vel and AoA peaks
for look_idx = 1:length(range_vals)
    find_idx(look_idx) = find(range_detect == range_vals(look_idx));
end
vel_vals = vel_detect(find_idx);
ang_vals = ang_detect(find_idx);

%display final results
fprintf('There are %d objects detected\n', length(range_vals));
for print_idx=1:length(range_vals)
    fprintf('Object %d has a range of %.2f m, a velocity of %.2f m/s, and an angle of attack of %.2f°.\n',print_idx, range_vals(print_idx),vel_vals(print_idx),ang_vals(print_idx)); 
end
