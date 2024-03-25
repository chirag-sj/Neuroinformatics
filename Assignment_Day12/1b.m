nyquist = EEG.srate/2;
lower_filter_bound = 4; % Hz
upper_filter_bound = 10; % Hz
transition_width   = 0.2;
filter_order       = round(3*(EEG.srate/lower_filter_bound));

% create the filter shape (this is explained more in the text around figure 14.4)
ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(filter_order,ffrequencies,idealresponse);

