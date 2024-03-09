
load('../sampleEEGdata (1).mat')
%%%%
% Parameters
frequencies = 2:5:30; % Frequencies from 2 Hz to 30 Hz in five steps
num_electrodes = size(EEG.data, 1); % Total number of electrodes

% Preallocate convolution results matrix
convolution_results = zeros(length(frequencies), EEG.pnts, num_electrodes);

% Loop through each frequency
for i = 1:length(frequencies)
    % Create Morlet wavelet for current frequency
    frequency = frequencies(i);
    time = -1:1/EEG.srate:1;
    s = (4/(2*pi*frequency))^2;
    wavelet = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*s)/frequency);
    
    % FFT parameters
    n_wavelet = length(wavelet);
    n_data = EEG.pnts;
    n_convolution = n_wavelet+n_data-1;
    half_of_wavelet_size = (length(wavelet)-1)/2;
    
    % FFT of wavelet
    fft_wavelet = fft(wavelet,n_convolution);
    
    % Loop through each electrode
    for electrode = 1:num_electrodes
        % FFT of EEG data for current electrode and trial 1
        fft_data = fft(squeeze(EEG.data(electrode,:,1)),n_convolution);
        
        % Convolution
        convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(s);
        
        % Cut off edges
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        
        % Store convolution result
        convolution_results(i, :, electrode) = convolution_result_fft;
    end
end

% Plot Morlet wavelets
figure;
for i = 1:length(frequencies)
    subplot(length(frequencies), 1, i);
    wavelet = exp(2*1i*pi*frequencies(i).*time) .* exp(-time.^2./(2*s)/frequencies(i));
    plot(time, real(wavelet)); % Plot the real part of the wavelet
    title(['Morlet Wavelet for Frequency ' num2str(frequencies(i)) ' Hz']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end
sgtitle('Morlet Wavelets for Different Frequencies');

% Plot convolution results for 1st trial and 47th electrode
figure;
for i = 1:length(frequencies)
    subplot(length(frequencies), 1, i);
    plot(real(convolution_results(i, :, 47))); % Plot convolution result for 47th electrode
    title(['Convolution Result for Frequency ' num2str(frequencies(i)) ' Hz']);
    xlabel('Time (samples)');
    ylabel('Amplitude');
    grid on;
end
sgtitle('Convolution Results for Different Frequencies (1st Trial, 47th Electrode)');
%%

power_phase_matrix = zeros(EEG.pnts, length(frequencies), num_electrodes, 2); % 2: 1 for power, 1 for phase

% Loop through each frequency
for i = 1:length(frequencies)
    % Loop through each electrode
    for electrode = 1:num_electrodes
        % Extract convolution result
        convolution_result = squeeze(convolution_results(i, :, electrode));
        
        % Extract power and phase
        power_phase_matrix(:, i, electrode, 1) = abs(convolution_result); % Power
        power_phase_matrix(:, i, electrode, 2) = angle(convolution_result); % Phase
    end
end
%%
tm_pnt = nearest(180*EEG.srate/1000);

figure
for i = 1:length(frequencies)
    subplot(3,2,i)
    plot_topography({EEG.chanlocs.labels},power_phase_matrix(tm_pnt,i,:,1));
    title(['Frequency = ' num2str(frequencies(i)) 'Hz'])
end
sgtitle('Topographical Distribution of Power at t=180s')
%%
tm_pnt = nearest(180*EEG.srate/1000);

figure
for i = 1:length(frequencies)
    subplot(3,2,i)
    plot_topography({EEG.chanlocs.labels},power_phase_matrix(tm_pnt,i,:,2));
    title(['Frequency = ' num2str(frequencies(i)) 'Hz'])
end
sgtitle('Topographical Distribution of Phase at t=180s')
%%
tm_pnt = nearest(360*EEG.srate/1000);

figure
for i = 1:length(frequencies)
    subplot(3,2,i)
    plot_topography({EEG.chanlocs.labels},power_phase_matrix(tm_pnt,i,:,1));
    title(['Frequency = ' num2str(frequencies(i)) 'Hz'])
end
sgtitle('Topographical Distribution of Power at t=360')
%%
tm_pnt = nearest(360*EEG.srate/1000);

figure
for i = 1:length(frequencies)
    subplot(3,2,i)
    plot_topography({EEG.chanlocs.labels},power_phase_matrix(tm_pnt,i,:,2));
    title(['Frequency = ' num2str(frequencies(i)) 'Hz'])
end
sgtitle('Topographical Distribution of Phase at t=360')


%%
%%%%
% Parameters
frequencies = [2 4 8 16 32 64]; % Frequencies from 2 Hz to 30 Hz in five steps
num_electrodes = size(EEG.data, 1); % Total number of electrodes

% Preallocate convolution results matrix
convolution_results = zeros(length(frequencies), EEG.pnts, num_electrodes);

% Loop through each frequency
for i = 1:length(frequencies)
    % Create Morlet wavelet for current frequency
    frequency = frequencies(i);
    time = -1:1/EEG.srate:1;
    s = (4/(2*pi*frequency))^2;
    wavelet = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*s)/frequency);
    
    % FFT parameters
    n_wavelet = length(wavelet);
    n_data = EEG.pnts;
    n_convolution = n_wavelet+n_data-1;
    half_of_wavelet_size = (length(wavelet)-1)/2;
    
    % FFT of wavelet
    fft_wavelet = fft(wavelet,n_convolution);
    
    % Loop through each electrode
    for electrode = 1:num_electrodes
        % FFT of EEG data for current electrode and trial 1
        fft_data = fft(squeeze(EEG.data(electrode,:,1)),n_convolution);
        
        % Convolution
        convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(s);
        
        % Cut off edges
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        
        % Store convolution result
        convolution_results(i, :, electrode) = convolution_result_fft;
    end
end

% Plot Morlet wavelets
figure;
for i = 1:length(frequencies)
    subplot(length(frequencies), 1, i);
    wavelet = exp(2*1i*pi*frequencies(i).*time) .* exp(-time.^2./(2*s)/frequencies(i));
    plot(time, real(wavelet)); % Plot the real part of the wavelet
    title(['Morlet Wavelet for Frequency ' num2str(frequencies(i)) ' Hz']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end
sgtitle('Morlet Wavelets for Different Frequencies');

% Plot convolution results for 1st trial and 47th electrode
figure;
for i = 1:length(frequencies)
    subplot(length(frequencies), 1, i);
    plot(real(convolution_results(i, :, 47))); % Plot convolution result for 47th electrode
    title(['Convolution Result for Frequency ' num2str(frequencies(i)) ' Hz']);
    xlabel('Time (samples)');
    ylabel('Amplitude');
    grid on;
end
sgtitle('Convolution Results for Different Frequencies (1st Trial, 47th Electrode)');
%%

power_phase_matrix = zeros(EEG.pnts, length(frequencies), num_electrodes, 2); % 2: 1 for power, 1 for phase

% Loop through each frequency
for i = 1:length(frequencies)
    % Loop through each electrode
    for electrode = 1:num_electrodes
        % Extract convolution result
        convolution_result = squeeze(convolution_results(i, :, electrode));
        
        % Extract power and phase
        power_phase_matrix(:, i, electrode, 1) = abs(convolution_result); % Power
        power_phase_matrix(:, i, electrode, 2) = angle(convolution_result); % Phase
    end
end
%%
tm_pnt = nearest(180*EEG.srate/1000);

figure
for i = 1:length(frequencies)
    subplot(3,2,i)
    plot_topography({EEG.chanlocs.labels},power_phase_matrix(tm_pnt,i,:,1));
    title(['Frequency = ' num2str(frequencies(i)) 'Hz'])
end
sgtitle('Topographical Distribution of Power at t=180s')
%%
tm_pnt = nearest(180*EEG.srate/1000);

figure
for i = 1:length(frequencies)
    subplot(3,2,i)
    plot_topography({EEG.chanlocs.labels},power_phase_matrix(tm_pnt,i,:,2));
    title(['Frequency = ' num2str(frequencies(i)) 'Hz'])
end
sgtitle('Topographical Distribution of Phase at t=180s')
%%
tm_pnt = nearest(360*EEG.srate/1000);

figure
for i = 1:length(frequencies)
    subplot(3,2,i)
    plot_topography({EEG.chanlocs.labels},power_phase_matrix(tm_pnt,i,:,1));
    title(['Frequency = ' num2str(frequencies(i)) 'Hz'])
end
sgtitle('Topographical Distribution of Power at t=360')
%%
tm_pnt = nearest(360*EEG.srate/1000);

figure
for i = 1:length(frequencies)
    subplot(3,2,i)
    plot_topography({EEG.chanlocs.labels},power_phase_matrix(tm_pnt,i,:,2));
    title(['Frequency = ' num2str(frequencies(i)) 'Hz'])
end
sgtitle('Topographical Distribution of Phase at t=360')
