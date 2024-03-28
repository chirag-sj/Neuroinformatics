% load("../sampleEEGdata (1).mat");
% nyquist = EEG.srate/2;
% lower_filter_bound = 4; % Hz
% upper_filter_bound = 10; % Hz
% transition_width   = 0.2;
% filter_order       = round(3*(EEG.srate/lower_filter_bound));
% 
% % create the filter shape (this is explained more in the text around figure 14.4)
% ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
% idealresponse = [ 0 0 1 1 0 0 ];
% filterweights = firls(filter_order,ffrequencies,idealresponse);
% 
% electrode_indx = [23 35 64];
% filtered_data = zeros(length(electrode_indx),EEG.pnts,EEG.trials);
% hilbert_data = zeros(size(filtered_data));
% phase_data = zeros(size(filtered_data));
% for i=1:length(electrode_indx)
%     filtered_data(i,:,:) = filtfilt(filterweights,1,double(squeeze(EEG.data(electrode_indx(i),:,:))));
%     hilbert_data(i,:,:) = hilbert(filtered_data(i,:,:));
% end
% 
% figure
% for i=1:length(electrode_indx)
%     subplot(3,1,i)
%     plot(EEG.times,abs(mean(exp(1i*squeeze(angle(hilbert_data(i,:,:)))),2)))
%     set(gca,'xlim',[-200 1000])
%     xlabel('Time (ms)')
%     ylabel('ITPC')
% end

load("../sampleEEGdata (1).mat");
nyquist = EEG.srate/2;
electrode_indx = [23 35 64];
frequencies = logspace(log10(4), log10(40), 20);
s = logspace(log10(3), log10(10), length(frequencies)) ./ (2 * pi * frequencies);

baselinetime = [ -300 -100 ];
baseidx=dsearchn(EEG.times',baselinetime(1)):dsearchn(EEG.times',baselinetime(2));

itpc_matrix_all = zeros(length(electrode_indx),length(frequencies), EEG.pnts);
power_matrix_all = zeros(length(electrode_indx),length(frequencies), EEG.pnts);
for elecInd=1:3
    itpc_matrix = zeros(length(frequencies), EEG.pnts);
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
    
        itpc_matrix(fi, :) = abs(mean(exp(1i * angle(hilbert_data)), 2));
        power_matrix(fi,:) = mean(abs(hilbert_data).^2,2);
        power_matrix(fi,:) = 10*log10(power_matrix(fi,:)./mean(power_matrix(fi,baseidx),2));
    end
    itpc_matrix_all(elecInd,:,:) = itpc_matrix(:, :);
    power_matrix_all(elecInd,:,:) = power_matrix(:,:);
end

baselinetime = [ -300 -100 ];
baseidx=dsearchn(EEG.times',baselinetime(1)):dsearchn(EEG.times',baselinetime(2));
for figInd=1:length(electrode_indx)
    figure;
    subplot(2,1,1)
    contourf(EEG.times,frequencies,squeeze(itpc_matrix_all(figInd,:,:)),40,'linecolor','none')
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