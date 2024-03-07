%% figure 1

% impulse function (all zeros; 1 in the middle)
impfun = zeros(1,600);
impfun(50)=1;
impfun(45:55)=1;

% kernel = [1 .8 .6 .4 .2]; 
% kernel = [.2 .9 1.0 .9 .2];
kernel = [1. .6 .4 .3 .25];


% matlab's convolution function
matlab_conv_result = conv(impfun,kernel,'same');

figure
% plot the signal (impulse or boxcar)
subplot(311)
plot(impfun)
set(gca,'ylim',[-.1 1.1])

% plot the kernel
subplot(312)
plot(kernel,'.-')
set(gca,'xlim',[0 100],'ylim',[-.1 1.1])

% plot the result of convolution
subplot(313)
plot(matlab_conv_result)
set(gca,'xlim',[0 100],'ylim',[-.1 3.6])

%% figure 2

% data that we'll use for convolution (must be zero-padded).
% dat4conv = [zeros(1,length(kernel)-1) impfun zeros(1,length(kernel)-1) ];
load_data = load("../sampleEEGdata (1).mat");
data_str = load_data.EEG;

% find(ismember(data_str.chanlocs.labels,"FCz"))
% find(data_str.chanlocs.labels=='FCz')
% strcmp(data_str.chanlocs.labels='FCz');

fcz_ind = find(strcmp({data_str.chanlocs.labels},'FCz'));

dat4conv_noPad = normalize(data_str.data(fcz_ind,1:50,1));
dat4conv = [zeros(1,length(kernel)-1) dat4conv_noPad zeros(1,length(kernel)-1) ];
% used for cutting the result of convolution
half_of_kernel_size = ceil((length(kernel)-1)/2);

% initialize convolution output
convolution_result = zeros(1,length(dat4conv_noPad)+length(kernel)-1);

% run convolution (note that kernel is flipped backwards)
for ti=1:length(convolution_result)-half_of_kernel_size
    convolution_result(ti) = sum(dat4conv(ti:ti+length(kernel)-1).*kernel(end:-1:1));
end

% cut off edges
convolution_result = convolution_result(half_of_kernel_size+1:end-half_of_kernel_size);
% Note: In the book figure there was a bug on the previous line ("+1" was typed 
%       as "-1"), which incorrectly showed the convolution as being a few points
%       too far to the right. Thanks to Thijs Perenboom for catching that!

figure
subplot(3,1,1)
plot(kernel)
title('Exponential Kernel')
subplot(3,1,2)
plot(dat4conv)
title('FCz Electrode Recording')
subplot(3,1,3)
plot(convolution_result)
title('Convolution Result')

%% Figure 3
kernel = [.2 .9 1.0 .9 .2];
dat4conv_noPad = normalize(data_str.data(fcz_ind,1:50,1));
dat4conv = [zeros(1,length(kernel)-1) dat4conv_noPad zeros(1,length(kernel)-1) ];
% used for cutting the result of convolution
half_of_kernel_size = ceil((length(kernel)-1)/2);

% initialize convolution output
convolution_result = zeros(1,length(dat4conv_noPad)+length(kernel)-1);

% run convolution (note that kernel is flipped backwards)
for ti=1:length(convolution_result)-half_of_kernel_size
    convolution_result(ti) = sum(dat4conv(ti:ti+length(kernel)-1).*kernel(end:-1:1));
end

% cut off edges
convolution_result = convolution_result(half_of_kernel_size+1:end-half_of_kernel_size);
% Note: In the book figure there was a bug on the previous line ("+1" was typed 
%       as "-1"), which incorrectly showed the convolution as being a few points
%       too far to the right. Thanks to Thijs Perenboom for catching that!

figure
subplot(3,1,1)
plot(kernel)
title('Inverted U Kernel')
subplot(3,1,2)
plot(dat4conv)
title('FCz Electrode Recording')
subplot(3,1,3)
plot(convolution_result)
title('Convolution Result')


%% 3. Based on visual inspection, what is the effect of convolving the EEG data with these two kernels?
% Both the filters are low pass filter which smoothen the EEG data as
% visible in the figure. Based on the figures we see that the gaussian
% filter is better at smoothing than its exponential counterpart. This is
% because exponential function leaves out some strong edges since nearer
% values are weighted more than farther ones thus distant ones have lower
% weightage.

%% Figure 4
kernel = [.2 .9 1.0 .9 .2];
dat4conv_noPad = normalize(data_str.data(fcz_ind,:,1));
% Define the signal
signal_length = length(dat4conv_noPad);
signal = dat4conv_noPad; % Example random signal
% Pad the kernel and the signal to the same length
padded_kernel = [kernel, zeros(1, signal_length - length(kernel))];
padded_signal = signal;%, zeros(1, length(kernel) - 1)];
% Perform FFT on the kernel and the signal
fft_kernel = fft(padded_kernel);
fft_signal = fft(padded_signal);
length(fft_kernel)
length(fft_signal)
% Perform the convolution by element-wise multiplication in the frequency domain
convolved = fft_kernel .* fft_signal;
% Inverse FFT to get the convolution result
result = ifft(convolved);
% Trim the result to the original signal length
result = result(1:signal_length);

figure
subplot(3,1,1)
plot(kernel)
title('Inverted U Kernel')
subplot(3,1,2)
plot(signal)
title('FCz Electrode Recording')
subplot(3,1,3)
plot(result)
title('Convolution Result using FFT')

%% Figure 5
kernel = [1. .6 .4 .3 .25];
dat4conv_noPad = normalize(data_str.data(fcz_ind,:,1));
% Define the signal
signal_length = length(dat4conv_noPad);
signal = dat4conv_noPad; % Example random signal
% Pad the kernel and the signal to the same length
padded_kernel = [kernel, zeros(1, signal_length - length(kernel))];
padded_signal = signal;%, zeros(1, length(kernel) - 1)];
% Perform FFT on the kernel and the signal
fft_kernel = fft(padded_kernel);
fft_signal = fft(padded_signal);
length(fft_kernel)
length(fft_signal)
% Perform the convolution by element-wise multiplication in the frequency domain
convolved = fft_kernel .* fft_signal;
% Inverse FFT to get the convolution result
result = ifft(convolved);
% Trim the result to the original signal length
result = result(1:signal_length);

figure
subplot(3,1,1)
plot(kernel)
title('Exponential Kernel')
subplot(3,1,2)
plot(signal)
title('FCz Electrode Recording')
subplot(3,1,3)
plot(result)
title('Convolution Result using FFT')