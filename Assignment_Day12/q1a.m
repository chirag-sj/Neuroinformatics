%% prereqs

%% Figure 18.1

figure

% 1/f function
c = 1;
x = 1;

plot(c./(1:100).^x)

%% Figure 18.2

% load sample EEG data
load("../sampleEEGdata (1).mat")

% wavelet parameters
min_freq = 2;
max_freq = 128;
num_frex = 30;

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet     = length(time);
n_data        = EEG.pnts;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));
wavelet_cycles= 4; 


electrode_indx = [23 35 64];
rand_tpnts = [10 50 150 300 640];
% FFT of data (note: this doesn't change on frequency iteration)
req_data = zeros(length(electrode_indx),length(rand_tpnts),EEG.trials);
for trial_ind=1:EEG.trials
    for elecInd=1:length(electrode_indx)
        hilbert_sig = hilbert(EEG.data(elecInd,:,trial_ind));
        inst_phase_all = angle(hilbert_sig);
        for tpnts=1:length(rand_tpnts)
            inst_phase_tpnts = inst_phase_all(tpnts);
            req_data(electrode_indx(elecInd),tpnts,trial_ind) = inst_phase_tpnts;
%             fft_data = EEG.data(electrode_indx(elecInd),rand_tpnts(tpnts),trial_ind);
%             disp(fft_data)
        end
    end
end
figure
elec_name_array = string({EEG.chanlocs.labels});
elecInd = electrode_indx(3); %select the electrode you want
% for figInd=1:15
%     subplot(3,5,figInd)
% for k=1:3

for j=1:5
    subplot(3,5,j)
    res_x = 0;
    res_y = 0;
    for i=1:EEG.trials
        % Compute x and y components of the vector
        x_component = cos(req_data(elecInd,j,i));   
        y_component = sin(req_data(elecInd,j,i));
        % Plot the unit vector starting from the origin
        quiver(0, 0, x_component, y_component, 'Color', 'b', 'LineWidth', 2);
        hold on;
%         res_x=res_x+x_component;
%         res_y=res_y+y_component;
%         if i==EEG.trials
%             quiver(0, 0, res_x, res_y, 'Color', 'r', 'LineWidth', 2);
%         end
    end
    axis equal;
    grid on;
    xlabel('X');
    ylabel('Y');
    title(['Timestamp:',num2str(rand_tpnts(j)),'Electrode',elec_name_array(elecInd)]);
end
% end
% end
