clc;
clear;
close all;

%% global parameters %%

%user interface
age_input = input('Age: ');
nervous = input('Filter type, 0 for strict, 1 for loose: ');
max_heart = 220-age_input;

%constants
N=1024; %FFT size
min_heart = 50; %minimum heart rate

switch nervous
    case 0
        heart_filter_cut = [0.83 1.67]; %Inital heart filter cutoff frequencies - > [50 100] breaths per min
    case 1
        heart_filter_cut = [0.83 2]; %Inital heart filter cutoff frequencies - > [50 120] breaths per min
end
    
breath_filter_cut = [0.2 0.33]; %Inital breath filter cutoff frequencies - > [12 20] breaths per min
filter_Fs = 20; %bandpass sampling frequency

%radar parameters (configurable) - provided by Dr. He
numADCSamples = 100;  % number of ADC samples per chirp
chirploop = 2; %number of chirps per frame 
frames = 400; %number of frames
Tr = 50e-6;  % Individual Chirp duration / Ramp end time
slope = 63e12; % Chirp slope 
Fs = 4e6; %Sampling Rate
Tf=0.05; % Frame Duration

%calculated values
numChirps = chirploop*frames; %number of total chirps
freq=(1:1:N/2-1)/Tf/N; %frequency spectrum for FFT plots

%file reading, raw data extraction, channel summation, and reshaping
fileName = 'C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc\adc_data.bin';
retVal = readDCA1000(fileName); %provided by Dr.He

%sum all channels to increase SNR
RXdata = zeros(numADCSamples,numChirps);
arr_dims = size(retVal);
for chan_num = 1:arr_dims(1)
    RXdata = RXdata + reshape(retVal(chan_num,:),numADCSamples,numChirps);
end

figure(1);
title('Radar Single Chirp Signal');
hold on
plot(db(real(RXdata(:,1)')),'DisplayName','Magnitude');
plot(db(imag(RXdata(:,1)')),'DisplayName','Phase');
xlabel('Magnitude (dB)');
ylabel('Sample #')
hold off
legend

%range fft
range_win = hamming(numADCSamples);
range_fft = fft(RXdata.*range_win);

figure(2);
plot(db(abs(range_fft)));
title('Range FFT Across All Chirps');
xlabel('Range (units)');
ylabel('Magnitude (dB)');

%angle fft in peak range bin
peak = zeros(1,frames);
for k = 1:frames
    for j = 1:numADCSamples
        if abs(range_fft(j, k)) == max(abs(range_fft(j,k)))
            peak(:, k) = range_fft(j, k);
        end
    end
end

%separate real and imaginary components of signal
data_real = zeros(1,frames);
data_imag = zeros(1,frames);
for k = 1:frames
    data_real(:, k) = real(peak(:, k));
    data_imag(:, k) = imag(peak(:, k));
end

%phase extraction and unwrap
signal_phase = zeros(1,frames);
for k = 1:frames
    signal_phase(:, k) = atan(data_imag(:, k) / data_real(:, k));
end

unwrap_phase = unwrap(signal_phase);

figure(3);
plot(unwrap_phase);
title('Unwrapped Signal Phase Extract');
xlabel('Sample #');
ylabel('Phase (°)');

%phase difference
phase_diff = diff(unwrap_phase);

figure(4);
plot(phase_diff);
title('Phase Difference');
xlabel('Sample #');
ylabel('Phase (°)');

%denoising
impulse_remove = medfilt1(phase_diff); %impulse noise filtering CHECK
%LATER TO SEE IF DELETE
ica_phase = fastICA(impulse_remove,1); %ICA separation - provided by Dr. He
init_filt = bandpass(ica_phase, [breath_filter_cut(1) heart_filter_cut(2)], filter_Fs); %initial filtering - gets rid of lowest and highest frequencies


figure(5);
plot(init_filt);
title('Denoised Signal');
xlabel('Sample #');
ylabel('Phase (°)');


%vitals separation ,fft, and normalization/scraping -> filter, fft + norm,
%account for symmetry, double energy
breath_filtered = bandpass(init_filt,[breath_filter_cut],filter_Fs); 
breath_fft = abs((fft(breath_filtered,N))/N); 
breath_half = breath_fft(1:N/2+1); 
breath_ener = 2*breath_half(2:end-1); 

heart_filtered = bandpass(init_filt,heart_filter_cut,filter_Fs); 
heart_fft = abs((fft(heart_filtered,N))/N); 
heart_half = heart_fft(1:N/2+1); 
heart_ener = 2*heart_half(2:end-1); 

figure(6);
plot(freq,db(breath_ener));
title('Respiratory Rate Extraction');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%cross reference with harmonic series - provided by Dr. He
%find peaks in heart-fft and 1st harmonic
[heart_peaks,heart_peaksnum]=findpeaks_custom(heart_ener,heart_filter_cut(1),heart_filter_cut(2),N,Tf); 
[heart_harmonic_peaks,heart_harmonic_peaksnum]=findpeaks_custom(heart_ener,heart_filter_cut(1)*2,heart_filter_cut(2)*2,N,Tf); 

%convert idx to freq
heart_peaks = heart_peaks/N/Tf;
heart_harmonic_peaks = heart_harmonic_peaks/N/Tf;

%dimensions for original data and harmonic
[heart_peaks_row,heart_peaks_column] = size(heart_peaks);
[heart_harmonic_peaks_row,heart_harmonic_peaks_column] = size(heart_harmonic_peaks);

%validation and strengthening
for i=1:heart_peaks_column
    if max(heart_ener)-heart_ener(round(heart_peaks(i)*N*Tf)+1)<0.3 
        for j=1:heart_harmonic_peaks_column
            if heart_harmonic_peaks(j)/heart_peaks(i)==2
                heart_ener(round(heart_peaks(i)*N*Tf)+1)=2*heart_ener(round(heart_peaks(i)*N*Tf)+1);
            end
        end
    end
end


figure(7);
plot(freq,db(heart_ener));
title('Harmonically Strengthened Heart Rate Extraction ');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%final predictions -> find peak mag, find related freq ind, convert to bpm
max_mag_breath = max(findpeaks(breath_ener));
max_ind_breath = find(breath_ener==max_mag_breath);
breath_predict = 60*freq(max_ind_breath);

max_mag_heart = max(findpeaks(heart_ener));
max_ind_heart = find(heart_ener==max_mag_heart);
heart_predict = 60*freq(max_ind_heart);

%output
fprintf('Respiratory Rate: %.2f Breaths Per Minute\n',breath_predict);
fprintf('Heart Rate: %.2f Beats Per Minute\n',heart_predict);

switch true
    case heart_predict > max_heart
        fprintf('Your heart rate is concerningly high, above your age-predicted maximum heart rate of %d BPM.\n',max_heart);
    case heart_predict < min_heart
        fprintf('Your heart rate is concerningly low, %d BPM below the average minimum resting heart rate of 50 .\n',min_heart-heart_predict);
end
