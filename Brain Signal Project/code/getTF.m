function [dbconverted] = getTF(trial_data)

% subject_ind=67;
% electrode='FCz';
% condition_val=1;
% [SZsub_data,demographic_sub,ERPdata_sub,mergedTrialData_sub,time_sub] = loadData(subject_ind);
time_sub = readtable('../archive/time.csv');
% trial_data = [];
% % hilbert_data = [];
% for i=1:100
%     if isempty(SZsub_data(SZsub_data.trial==i & SZsub_data.condition==condition_val,:).(electrode))
%         continue
%     elseif isempty(trial_data)
%         trial_data = detrend(SZsub_data(SZsub_data.trial==i & SZsub_data.condition==condition_val,:).(electrode)');
%     else
%         sub_trial_data = detrend(SZsub_data(SZsub_data.trial==i & SZsub_data.condition==condition_val,:).(electrode)');
%         trial_data = cat(1, trial_data,sub_trial_data);
%     end
% end
min_freq = 2;
max_freq = 64;
num_frex = 30;

% other wavelet parameters
srate=1024;
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1.5:1/srate:1.5;
half_of_wavelet_size = (length(time)-1)/2;

pnts = 3072;
n_wavelet     = length(time);
n_data        = pnts;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));
wavelet_cycles= 4; 

% FFT of data (note: this doesn't change on frequency iteration)
fft_data = fft(squeeze(trial_data(1,:)),n_conv_pow2);

% initialize output time-frequency data
tf_data = zeros(length(frequencies),pnts);

for fi=1:length(frequencies)
    
    % create wavelet and get its FFT
    wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
    fft_wavelet = fft(wavelet,n_conv_pow2);
    
    % run convolution
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
    convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % put power data into time-frequency matrix
    tf_data(fi,:) = abs(convolution_result_fft).^2;
end

% plot results
% ytickskip = 2:4:num_frex; % This will be explained in the text.
% figure
% imagesc(time_sub.time_ms,[],tf_data)
% set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','clim',[-12 12]) %'xlim',[-500 1500],
% xlabel('Time (ms)'), ylabel('Frequency (Hz)')


% define baseline period
baselinetime = [ -500 -200 ]; % in ms

% convert baseline window time to indices
[~,baselineidx(1)]=min(abs(time_sub.time_ms-baselinetime(1)));
[~,baselineidx(2)]=min(abs(time_sub.time_ms-baselinetime(2)));


baseline_power = mean(tf_data(:,baselineidx(1):baselineidx(2)),2);
dbconverted = 10*log10( bsxfun(@rdivide,tf_data,baseline_power));

% figure
% contourf(time_sub.time_ms,frequencies,dbconverted,40,'linecolor','none')
% set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','clim',[-12 12])