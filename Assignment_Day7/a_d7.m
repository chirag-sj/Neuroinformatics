% Create a family of Morlet wavelets ranging in frequency from 2 Hz to 30 Hz.
% [lb,ub,n,fb,fc] = -1,1,20,1,2;
[psi_2,x] = cmorwavf(-1,1,30,1,2);
[psi_4,~] = cmorwavf(-1,1,30,1,4);
[psi_8,~] = cmorwavf(-1,1,30,1,8);
[psi_15,~] = cmorwavf(-1,1,30,1,15);
[psi_30,~] = cmorwavf(-1,1,30,1,30);
figure
subplot(5,1,1)
plot(x,real(psi_2))
subplot(5,1,2)
plot(x,real(psi_4))
subplot(5,1,3)
plot(x,real(psi_8))
subplot(5,1,4)
plot(x,real(psi_15))
subplot(5,1,5)
plot(x,real(psi_30))
 
% Convolve each Morlet wavelet with EEG data from all trials from the electrode we used in the previous day. 
% Apply the Matlab function real to the convolution result, as in convol_result=real(convol_result). 
% This will return the EEG data bandpass filtered at the peak frequency of the wavelet.
eegdata = load("../sampleEEGdata (1).mat");
data = normalize(eegdata.EEG.data,3);
conv_data_2 = zeros(size(data));
conv_data_4 = zeros(size(data));
conv_data_8 = zeros(size(data));
conv_data_15 = zeros(size(data));
conv_data_30 = zeros(size(data));

for elec = 1:64
    for trl = 1:99
        conv_data_2(elec,:,trl) = conv(data(elec,:,trl),psi_2,'same');
        conv_data_4(elec,:,trl) = conv(data(elec,:,trl),psi_4,'same');
        conv_data_8(elec,:,trl) = conv(data(elec,:,trl),psi_8,'same');
        conv_data_15(elec,:,trl) = conv(data(elec,:,trl),psi_15,'same');
        conv_data_30(elec,:,trl) = conv(data(elec,:,trl),psi_30,'same');
    end
end
figure
plot(real(data(47,:,1)))

figure
subplot(5,1,1)
plot(real(conv_data_2(47,:,1)))
title('FCz electrode 2Hz Trial 1')
subplot(5,1,2)
plot(real(conv_data_4(47,:,1)))
title('FCz electrode 4Hz Trial 1')
subplot(5,1,3)
plot(real(conv_data_8(47,:,1)))
title('FCz electrode 8Hz Trial 1')
subplot(5,1,4)
plot(real(conv_data_15(47,:,1)))
title('FCz electrode 15Hz Trial 1')
subplot(5,1,5)
plot(real(conv_data_30(47,:,1)))
title('FCz electrode 30Hz Trial 1')

% 
% Average the result of convolution over all trials and plot an ERP corresponding to each wavelet frequency. 
% Each frequency should be in its own subplot. Plot the broadband ERP (without any convolution) too.
avg_2 = mean(conv_data_2,3);
avg_4 = mean(conv_data_4,3);
avg_8 = mean(conv_data_8,3);
avg_15 = mean(conv_data_15,3);
avg_30 = mean(conv_data_30,3);
% size(avg_2)
figure
subplot(5,1,1)
plot(real(avg_2(47,:)'))
title('FCz ERP 2Hz')
subplot(5,1,2)
plot(real(avg_4(47,:)'))
title('FCz ERP 4Hz')
subplot(5,1,3)
plot(real(avg_8(47,:)'))
title('FCz ERP 8Hz')
subplot(5,1,4)
plot(real(avg_15(47,:)'))
title('FCz ERP 15Hz')
subplot(5,1,5)
plot(real(avg_30(47,:)'))
title('FCz ERP 30Hz')


% How do the wavelet-convolved ERPs compare with the broadband ERP?
gauss_data = zeros(size(data));
sigma = 1;
u=0;
t1 = linspace(-1,1,30);
Exp_comp = -((t1-u).^2)/(2*sigma*sigma);%tn/(sigma*sqrt(2))
G = exp(Exp_comp)/((sqrt(2*pi))*sigma);
for elec = 1:64
    for trl = 1:99
        gauss_data(elec,:,trl) = conv(data(elec,:,trl),G,'same');
    end
end

broadERP = mean(data,3);
gaussBroudERP = mean(gauss_data,3);
figure
subplot(2,1,1)
plot(broadERP(47,:))
title('Without Gaussian Filtering')
subplot(2,1,2)
plot(gaussBroudERP(47,:))
title('With Gaussian Filtering')

%The broadband ERP is not as interpretable as the wavelet based ERPs. The
%wavelet based ERPs amplifies the pattern which are specific to those
%frequency ranges.
