%% Figure 19.12
load("../sampleEEGdata (1).mat");
centerfreq  = 6;
channel2use = 'po7';
times2save  = -200:50:1200;

% initialize matrix to store RTs
rts = zeros(size(EEG.epoch));

for ei=1:EEG.trials
    
    % find which event is time=0, and take the latency of the event thereafter.
    time0event = find(cell2mat(EEG.epoch(ei).eventlatency)==0);
    
    % use try-catch in case of no response
    try
        rts(ei) = EEG.epoch(ei).eventlatency{time0event+1};
    catch me;
        rts{ei} = NaN;
    end
end


% % definte convolution parameters
% time          = -EEG.pnts/EEG.srate/2:1/EEG.srate:EEG.pnts/EEG.srate/2-1/EEG.srate;
% n_wavelet     = EEG.pnts;
% n_data        = EEG.pnts*EEG.trials;
% n_convolution = n_wavelet+n_data-1;
% n_conv_pow2   = pow2(nextpow2(n_convolution));
% 
% % get FFT of data and wavelet
% eegfft  = fft(reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,n_data),n_conv_pow2);
% wavefft = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*((4/(2*pi*centerfreq))^2))),n_conv_pow2);
% 
% % convolution
% eegconv = ifft(wavefft.*eegfft);
% eegconv = eegconv(1:n_convolution);
% eegconv = reshape(eegconv(floor((EEG.pnts-1)/2):end-1-ceil((EEG.pnts-1)/2)),EEG.pnts,EEG.trials);
% 
% phase_angles = angle(eegconv);
% 
% % initialize
% itpc    = zeros(size(times2save));
% witpc   = zeros(size(times2save));
% witpc_z = zeros(size(times2save));
% 
% for ti=1:length(times2save)
%     
%     % find index for this time point
%     [junk,timeidx] = min(abs(EEG.times-times2save(ti)));
%     
%     % ITPC is unmodulated phase clustering
%     itpc(ti) = abs(mean(exp(1i*phase_angles(timeidx,:))));
%     
%     % wITPC is rts modulating the length of phase angles
%     witpc(ti) = abs(mean(rts.*exp(1i*phase_angles(timeidx,:))));
%     
%     % permutation testing
%     perm_witpc = zeros(1,1000);
%     for i=1:1000
%         perm_witpc(i) = abs(mean(rts(randperm(EEG.trials)).*exp(1i*phase_angles(timeidx,:))));
%     end
%     
%     witpc_z(ti) = (witpc(ti)-mean(perm_witpc))./std(perm_witpc);
% end
% 
% figure
% 
% % plot ITPC
% subplot(311)
% h=plotyy(EEG.times,abs(mean(exp(1i*phase_angles),2)),times2save,witpc);
% set(h,'xlim',[times2save(1) times2save(end)]);
% title([ 'ITPC and wITPC at ' channel2use ])
% legend({'ITPC';'wITPC'})
% 
% % plot wITPCz
% subplot(312)
% plot(times2save,witpc_z)
% hold on
% plot(get(gca,'xlim'),[0 0],'k')
% set(gca,'xlim',[times2save(1) times2save(end)])
% xlabel('Time (ms)'), ylabel('wITPCz')
% title([ 'wITPCz at ' channel2use ])
% 
% 
% subplot(325)
% plot(itpc,witpc,'.')
% xlabel('ITPC'), ylabel('wITPC')
% axis square
% 
% subplot(326)
% plot(itpc,witpc_z,'.')
% xlabel('ITPC'), ylabel('wITPCz')
% axis square


load("../sampleEEGdata (1).mat");
nyquist = EEG.srate/2;
electrode_indx = [23 35 64];
frequencies = logspace(log10(4), log10(40), 20);
s = logspace(log10(3), log10(10), length(frequencies)) ./ (2 * pi * frequencies);

baselinetime = [ -300 -100 ];
baseidx=dsearchn(EEG.times',baselinetime(1)):dsearchn(EEG.times',baselinetime(2));

witpc_matrix_all = zeros(length(electrode_indx),length(frequencies), EEG.pnts);
power_matrix_all = zeros(length(electrode_indx),length(frequencies), EEG.pnts);
for elecInd=1:3
    witpc_matrix = zeros(length(frequencies), EEG.pnts);
    power_matrix = zeros(length(frequencies), EEG.pnts);
    for fi = 1:length(frequencies)
        % Create filter
        lower_filter_bound = frequencies(fi) - 0.5*s(fi);
        upper_filter_bound = frequencies(fi) + 0.5*s(fi);
        transition_width = 0.2;
        filter_order = round(3 * (EEG.srate / lower_filter_bound));
        ffrequencies = [0 (1 - transition_width) * lower_filter_bound lower_filter_bound upper_filter_bound (1 + transition_width) * upper_filter_bound nyquist] / nyquist;
        idealresponse = [0 0 1 1 0 0];
        filterweights = firls(filter_order, ffrequencies, idealresponse);
    
        filtered_data = filtfilt(filterweights, 1, double(squeeze(EEG.data(electrode_indx(elecInd), :, :)))); % Assuming you're interested in the first electrode
    
        hilbert_data = hilbert(filtered_data);
    
        witpc_matrix(fi, :) = abs(mean(repmat(rts,640,1).*exp(1i * angle(hilbert_data)), 2));%rts.*
        power_matrix(fi,:) = mean(abs(hilbert_data).^2,2);
        power_matrix(fi,:) = 10*log10(power_matrix(fi,:)./mean(power_matrix(fi,baseidx),2));
    end
    witpc_matrix_all(elecInd,:,:) = witpc_matrix(:, :);
    power_matrix_all(elecInd,:,:) = power_matrix(:,:);
end

baselinetime = [ -300 -100 ];
baseidx=dsearchn(EEG.times',baselinetime(1)):dsearchn(EEG.times',baselinetime(2));
for figInd=1:length(electrode_indx)
    figure;
    subplot(2,1,1)
    contourf(EEG.times,frequencies,squeeze(witpc_matrix_all(figInd,:,:)),40,'linecolor','none')
    set(gca,'clim',[.1 .5],'xlim',[-200 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title('ITPC')
    
    subplot(2,1,2)
    contourf(EEG.times,frequencies,squeeze(power_matrix_all(figInd,:,:)),40,'linecolor','none')
    set(gca,'clim',[-1 1],'xlim',[-200 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title('DB-power ')
end