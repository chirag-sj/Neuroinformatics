% Compute the short-time FFT at each electrode and make topographical maps of thetaband (around 6 Hz) power and 
% alpha-band (around 10 Hz) power at 150 ms and 700 ms.
% create han window. apply fft
%% Figure 15.1
load('../sampleEEGdata (1).mat');
% timewin = 500; % in ms
% % convert ms to idx
% timewinidx = round(timewin/(1000/EEG.srate));
% % create hann taper function
% hann_win = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));
% % detrend data (useful to attentuate super-low frequency artifacts in FFT
% % from sampled data)
% d = detrend(EEG.data(:,:,16)'); %works on each column separately. 64x640 size matrix. So take transpose.
% % disp(size(d));
% figure
% subplot(311)
% plot(EEG.times,d(:,20))
% title('One trial of data') 
% % % comparing detrend with original
% % subplot(312)
% % plot(EEG.times,EEG.data(:,:,16))
% %the signal is to be analyzed at time = 150ms, 700ms
% [~,stime] = min(abs(EEG.times--150));
% subplot(323)
% plot(EEG.times(stime:stime+timewinidx-1),d(stime:stime+timewinidx-1,20),'b')
% hold on
% d_place = bsxfun(@times,d(stime:stime+timewinidx-1,:),hann_win');
% % size(bsxfun(@times,d(stime:stime+timewinidx-1,:),hann_win')) %bsxfun
% % function enables us to multiply (@times) timeseries corresponding to each
% % electrode with the Hann window.
% plot(EEG.times(stime:stime+timewinidx-1),d_place(:,20),'r')
% set(gca,'xlim',[-50 -50+timewin])
% title('One short-time window of data, windowed')
% % 
% dfft = fft(d_place); %size(d_place) is 128x64. Each column is ts corresponding to one electrode. FFT operates on columns.
% f    = linspace(0,EEG.srate/2,floor(length(hann_win)/2)+1); % frequencies of FFT
% subplot(313)
% plot(f(2:end),abs(dfft(2:floor(length(hann_win)/2)+1,20)).^2,'.-');
% title('power spectrum from that time window')
% set(gca,'xlim',[1 128],'ylim',[-1000 25000],'xtick',0:10:EEG.srate/2)
% ts
% % create TF matrix and input column of data at selected time point
% tf = zeros(floor(length(hann_win)/2),EEG.pnts);
% % tf(:,stime+timewinidx/2-10:stime+timewinidx/2+10,:) = repmat(abs(dfft(2:floor(length(hann_win)/2)+1,:))*2,1,1,21);
% tf(:,stime+timewinidx/2-10:stime+timewinidx/2+10) = repmat(abs(dfft(2:floor(length(hann_win)/2)+1,20))*2,1,21);
% figure
% imagesc(EEG.times,f,log10(tf(:,:)+1))
% set(gca,'clim',[-1 1]*4)
% 
% figure 
% plot_topography({EEG.chanlocs.labels},tf(1,1,:));

%% Q2
% Select one electrode and one frequency and compute power over time at that electrode and that frequency for all trials using 
elec_ind = 47;
center_freq = 10;
%% (1) complex wavelet convolution, 

frequency = center_freq;
time = -1:1/EEG.srate:1;
s = (4/(2*pi*frequency))^2;
wavelet = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*s)/frequency);
data_elec = EEG.data(elec_ind,:,:); %each column is 1 trial

% FFT parameters
n_wavelet = length(wavelet);
n_data = EEG.pnts;
n_convolution = n_wavelet+n_data-1;
half_of_wavelet_size = (length(wavelet)-1)/2;

% FFT of wavelet
fft_wavelet = fft(wavelet,n_convolution);
convolution_results = zeros(n_data,EEG.trials);
for trial_ind = 1:EEG.trials
    % FFT of EEG data for current electrode and trial 1
    fft_data = fft(data_elec(1,:,trial_ind),n_convolution);
    
    % Convolution
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(s);
    
    % Cut off edges
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Store convolution result
    convolution_results(:, trial_ind) = convolution_result_fft;
end

power_phase_matrix = zeros(EEG.pnts,EEG.trials, 2); % 2: 1 for power, 1 for phase
% Loop through each frequency
for i = 1:EEG.trials
    % Extract power and phase
    power_phase_matrix(:,i, 1) = power(abs(convolution_results(:,i)),2); % Power
    power_phase_matrix(:,i, 2) = angle(convolution_results(:,i)); % Phase
end

figure 
% for i=1:1
%     plot(EEG.times,power_phase_matrix(:,i,1),DisplayName=strcat('Trial',num2str(i)))
%     hold on
% end
pow_phase_mat = mean(power_phase_matrix,2);
plot(EEG.times,pow_phase_mat(:,1,1))
title('CMW')

%% (2) filter-Hilbert, 
%make the signal narrow band using filter then take hilbert transform, then
%get analytic signal and finally calculate power and phase.

center_freq   = 10;
freqspread    = 2; % Hz +/- the center frequency
transwid      = .10;
nyquist = EEG.srate/2;
electrode_ind = 47;
idealresponse = [0 0 1 1 0 0];
ffrequencies  = [ 0 (1-transwid)*(center_freq-freqspread) (center_freq-freqspread) (center_freq+freqspread) (1+transwid)*(center_freq+freqspread) nyquist ]/nyquist;

data2filter   = squeeze(double(EEG.data(electrode_ind,:,1)));
filterweights = firls(200,ffrequencies,idealresponse); % recompute without z-scoring

convol_result = conv(data2filter,filterweights,'same'); % could also use ifft(fft(data2filter...
analytic_signal = hilbert(convol_result);
power_of_analytic_sig = abs(analytic_signal);

figure
plot(EEG.times,power_of_analytic_sig)
xlabel('Time (ms)'), ylabel('Power (\muV^2)');
legend({'Power'})
title('filter-Hilbert')
%% (3) the short-time FFT, and 

timewin        = 400; % in ms, for stFFT
times2save     = -300:50:1000; % in ms
% channel2plot   = 'p7';
timepoint2plot = 200; % ms

% convert from ms to index
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end
timewinidx = round(timewin/(1000/EEG.srate));
% chan2useidx = strcmpi(channel2plot,{EEG.chanlocs.labels});
chan2useidx = 47;

% create hann taper
hann_win = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));

% define frequencies
frex = linspace(0,EEG.srate/2,floor(timewinidx/2)+1);

% initialize power output matrix
tf = zeros(length(frex),length(times2save));

% loop over time points and perform FFT
for timepointi=1:length(times2save)
    
    % extract time series data for this center time point
    tempdat = squeeze(EEG.data(chan2useidx,times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2),:)); % note: the 'mod' function here corrects for even or odd number of points
    
    % taper data (using bsxfun instead of repmat... note sizes of tempdat
    % and hann_win)
    taperdat = bsxfun(@times,tempdat,hann_win');
    
    fdat = fft(taperdat,[],1)/timewinidx; % 3rd input is to make sure fft is over time
    tf(:,timepointi) = mean(abs(fdat(1:floor(timewinidx/2)+1,:)).^2,2); % average over trials
end

chan2use = 'FCz';
frequency2plot = 10;  % in Hz

[junk,freq2plotidx]=min(abs(frex-frequency2plot));

% initialize ITPC output matrix
itpc = zeros(length(frex),length(times2save));

% loop over time points and perform FFT
for timepointi=1:length(times2save)
    
    % extract time series data for this center time point
    % (yes, yes, it's a long line of code. Perhaps you can understand it
    % better by breaking it up into several lines to separately identify
    % channel index and time window?)
    tempdat = squeeze(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2),:));
    
    % taper data (Note the difference between using repmat here and bsxfun
    % earlier in the code.)
    taperdat = tempdat.*repmat(hann_win',1,EEG.trials);
    
    fdat = fft(taperdat,[],1)/timewinidx; % 3rd input is to make sure fft is over time
    itpc(:,timepointi) = abs(mean(abs(fdat(1:floor(timewinidx/2)+1,:)).^2,2)); % average over trials
end


figure
plot(times2save,mean(itpc(freq2plotidx-2:freq2plotidx+2,:),1))
title('Short FFT')

% (4) multi-taper method.


channel2plot    = 'FCz';
frequency2plot  = 10;  % in Hz
timepoint2plot  = 200; % ms

nw_product      = 3;  % determines the frequency smoothing, given a specified time window
times2save      = -300:50:1000;
baseline_range  = [-200 -00];
timewin         = 400; % in ms

% convert time points to indices
times2saveidx = dsearchn(EEG.times',times2save'); 
timewinidx    = round(timewin/(1000/EEG.srate));

% find baselinetimepoints
baseidx = zeros(size(baseline_range));
[~,baseidx(1)] = min(abs(times2save-baseline_range(1)));
[~,baseidx(2)] = min(abs(times2save-baseline_range(2)));
% note that the following line is equivalent to the previous three
%baseidx = dsearchn(times2save',baseline_range');

% define tapers
tapers = dpss(timewinidx,nw_product); % note that in practice, you'll want to set the temporal resolution to be a function of frequency
% define frequencies for FFT
f = linspace(0,EEG.srate/2,floor(timewinidx/2)+1);

% find logical channel index
chanidx = strcmpi(channel2plot,{EEG.chanlocs.labels});

% initialize output matrix
multitaper_tf = zeros(floor(timewinidx/2)+1,length(times2save));

% loop through time bins
for ti=1:length(times2saveidx)
    
    % initialize power vector (over tapers)
    taperpow = zeros(floor(timewinidx/2)+1,1);
    
    % loop through tapers
    for tapi = 1:size(tapers,2)-1
        
        % window and taper data, and get power spectrum
        data      = bsxfun(@times,squeeze(EEG.data(chanidx,times2saveidx(ti)-floor(timewinidx/2)+1:times2saveidx(ti)+ceil(timewinidx/2),:)),tapers(:,tapi));
        pow       = fft(data,timewinidx)/timewinidx;
        pow       = pow(1:floor(timewinidx/2)+1,:);
        taperpow  = taperpow + mean(pow.*conj(pow),2);
    end
    
    % finally, get power from closest frequency
    multitaper_tf(:,ti) = taperpow/tapi;
end

% db-correct
db_multitaper_tf = 10*log10( multitaper_tf ./ repmat(mean(multitaper_tf(:,baseidx(1):baseidx(2)),2),1,length(times2save)) );


% plot time courses at one frequency band
figure
% subplot(121)
[junk,freq2plotidx]=min(abs(f-frequency2plot)); % can replace "junk" with "~"
plot(times2save,mean(log10(multitaper_tf(freq2plotidx-2:freq2plotidx+2,:)),1))
title('Multitaper')
% axis square
% set(gca,'xlim',[times2save(1) times2save(end)])

% subplot(122)
% [junk,time2plotidx]=min(abs(times2save-timepoint2plot));
% plot(f,log10(multitaper_tf(:,time2plotidx)))
% title([ 'Sensor ' channel2plot ', ' num2str(timepoint2plot) ' ms' ])

% Plot the trial-averaged results of these four time-frequency decomposition methods in different subplots of one figure.
% 
% Note that the scaling might be different because no baseline normalization has been applied.
% 
% Can you think of an appropriate baseline normalization/scaling?
% 
% How visually similar are the results from these three methods? If the results from the three methods are different, how are they different, and what parameters do you think you could change in the three methods to make the results look more or less similar?

