% subject_ind=1;
% electrode='FCz';
% condition_val=1;
% [sub_data,demographic_sub,ERPdata_sub,mergedTrialData_sub,time_sub] = loadData(subject_ind);
% trial_data = [];
% % hilbert_data = [];
% for i=1:100
%     if isempty(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).(electrode))
%         continue
%     elseif isempty(trial_data)
%         trial_data = detrend(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).(electrode)');
% %         hilbert_data = hilbert_data(trial_data);
%     else
%         sub_trial_data = detrend(sub_data(sub_data.trial==i & sub_data.condition==condition_val,:).(electrode)');
%         trial_data = cat(1, trial_data,sub_trial_data);
% %         hilbert_data = cat(1,hilbert_data,hilbert(sub_trial_data));
%     end
% end
% 
% hilbert_data = hilbert(trial_data');
% % plot(angle(hilbert_data(:,1)),Color=[1 0 0])
% % hold on;
% % plot(angle(hilbert(trial_data(1,:))),Color=[0 0 1])
% 
% % TF plot of ITPC
% n_wavelet     = 3072;
% n_data        = 3072*100;
% n_convolution = n_wavelet+n_data-1;
% n_conv_pow2   = pow2(nextpow2(n_convolution));
% 
% frequencies = logspace(log10(4),log10(40),20);
% s = logspace(log10(3),log10(10),length(frequencies))./(2*pi*frequencies);
% itpc = zeros(length(frequencies),3072);
% time    = -3072/1024/2:1/1024:3072/1024/2-1/1024;
% eegfft = fft(reshape(trial_data,1,[]),n_conv_pow2);
% for fi=1:length(frequencies)
%     % create wavelet
%     wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*(s(fi)^2)))/frequencies(fi);
%     % convolution
%     eegconv = ifft(fft(wavelet,n_conv_pow2).*eegfft);
%     eegconv = eegconv(1:n_convolution);
%     eegconv = reshape(eegconv(floor((3072-1)/2):end-1-ceil((3072-1)/2)),3072,100);
%     
%     % extract ITPC
%     itpc(fi,:) = abs(mean(exp(1i*angle(eegconv)),2));
% end

figure
contourf(time_sub.time_ms,frequencies,itpc,40,'linecolor','none')
set(gca,'clim',[-0.1 0.1])%,'xlim',[-200 1000]
xlabel('Time (ms)'), ylabel('Frequencies (Hz)')
