function [overthreshold_interval_first,ups,durations,first_peak,peak,down_crossing_index ,up_crossing_index, latency_threshold, thr_orig, integ_baseline,aver_fr,ejex] = get_latency_newest_ons(whole,base_spikes,num_base,alpha_gauss,sr,half_ancho_gauss,tmin_base,tmax_base,t_offset,tmin_tot,tmax_tot,num_ok,nstd,nstd_fact,over_threshold_time,t_down,lat_scale,thr_min)

wind_type = 'gaussian';
% wind_type = 'lapla';

% nstd = 5;   % number of standard deviations that the fire has to be over
% the baseline level of the neuron to be considered as the
% response

% firing rate calculation parameters
% sigma_gauss = 50;   % the width of the gaussian on the convolution window in ms
% t_offset = tmin_post;
% t_offset = 0;

sample_period = 1000/sr; % sample period for the spike list - window convolution in ms/sample
decimate_factor = 1; % decimation factor used after windowing the spikes times.
% alpha_gauss = 3; % ALPHA parameter for the gausswin function
% 2.5 => last window value ~ 5%
% 3   => last window value ~ 1%
N_gauss = 2*round(half_ancho_gauss/sample_period)+1; % Number of points of the gaussian window
% Gauss_wind_time = N_gauss * sample_period; % time extension of the gaussian window in ms.

ejex = tmin_tot-half_ancho_gauss:sample_period:tmax_tot+half_ancho_gauss;  %should be the same length as spike_timeline

% creation of the time window for integrate spikes
if strcmp(wind_type,'gaussian')
    int_window = gausswin(N_gauss, alpha_gauss);
elseif strcmp(wind_type,'lapla')
    tempos = ((1:N_gauss)-ceil(N_gauss/2))/floor(N_gauss/2)*half_ancho_gauss;
    b_lapla = half_ancho_gauss/alpha_gauss;
    int_window = exp(-abs(tempos)/b_lapla)/2/b_lapla;
end
if lat_scale==1
    int_window = 1000*int_window/sum(int_window)/sample_period;
end

times_ind_base = find(ejex>tmin_base,1)-1:find(ejex>tmax_base,1);
times_ind_post = find(ejex>t_offset,1)-1:length(ejex);

total_spikes = sort(whole);
spike_timeline = hist(total_spikes,(tmin_tot-half_ancho_gauss:sample_period:tmax_tot+half_ancho_gauss))/num_ok;
n_spike_timeline = length(spike_timeline); %should be the same length as ejex
integ_timeline_stim = conv(spike_timeline, int_window);
integ_timeline_stim_cut = integ_timeline_stim(round(half_ancho_gauss/sample_period)+1:n_spike_timeline+round(half_ancho_gauss/sample_period));
total_base_spikes = sort(base_spikes);
base_spike_timeline = hist(total_base_spikes,(tmin_tot-half_ancho_gauss:sample_period:tmax_tot+half_ancho_gauss))/num_base;
n_base_spike_timeline = length(base_spike_timeline); %should be the same length as ejex
integ_timeline_base = conv(base_spike_timeline, int_window);

integ_timeline_base_cut = integ_timeline_base(round(half_ancho_gauss/sample_period)+1:n_base_spike_timeline+round(half_ancho_gauss/sample_period));
integ_baseline = mean(integ_timeline_base_cut(times_ind_base));
integ_std_dev = std(integ_timeline_base_cut(times_ind_base));
thr_orig = (nstd*integ_std_dev)+integ_baseline;
latency_threshold = max(thr_orig,min(thr_min,nstd_fact*thr_orig)); % fixed to correct bug with zero baseline.
timeline = integ_timeline_stim_cut(times_ind_post); % starts looking after stimulus onset

over_threshold = timeline > latency_threshold;
over_th_shift = circshift(over_threshold,[0 1]); % the first component should be forced to 0 to avoid problems if the last one is 1
over_th_shift(1)=0;
over_threshold(1)=0;
crossing_points = over_threshold - over_th_shift;
up_crossing_index = find(crossing_points == 1);
down_crossing_index = find(crossing_points == -1);
overthreshold_interval_index = [];
if ~isempty(up_crossing_index)
    if isempty(down_crossing_index)
        down_crossing_index = length(crossing_points);
    end
    if down_crossing_index(1)<up_crossing_index(1)
        down_crossing_index = down_crossing_index(2:end);
    end
    if length(up_crossing_index) == length(down_crossing_index) + 1    % if length(up_crossing_index) = length(down_crossing_index) + 1, it means that the last interval is length(crossing_points) + 1 - up_crossing_index(end)
        if length(up_crossing_index) > 1
            overthreshold_interval_index = down_crossing_index - up_crossing_index(1:end-1);
        end
        overthreshold_interval_index = [overthreshold_interval_index,length(crossing_points)-up_crossing_index(end)];
        down_crossing_index = [down_crossing_index length(crossing_points)];
        if length(up_crossing_index) > 1
            saca = (up_crossing_index(2:end)-down_crossing_index(1:end-1))*sample_period<t_down;
            down_crossing_index(find(saca))=[];
            up_crossing_index(find(saca)+1)=[];
            overthreshold_interval_index = down_crossing_index - up_crossing_index;
        end
        overthreshold_interval_first = find (overthreshold_interval_index >= (over_threshold_time/sample_period),1);
    else
        overthreshold_interval_index = down_crossing_index - up_crossing_index;
        if length(up_crossing_index) > 1
            saca = (up_crossing_index(2:end)-down_crossing_index(1:end-1))*sample_period<t_down;
            down_crossing_index(find(saca))=[];
            up_crossing_index(find(saca)+1)=[];
            overthreshold_interval_index = down_crossing_index - up_crossing_index;
        end
        overthreshold_interval_first = find (overthreshold_interval_index >= (over_threshold_time/sample_period),1);
    end
else
    overthreshold_interval_first = [];
end

if ~isempty(overthreshold_interval_index)
    ups = ejex(times_ind_post(up_crossing_index));
    durations = overthreshold_interval_index*sample_period;
else
    ups = -1;
    durations = -1;
end

first_peak = -1;

if ~isempty(overthreshold_interval_first)
    inds_dura = times_ind_post(up_crossing_index(overthreshold_interval_first)):times_ind_post(down_crossing_index(overthreshold_interval_first));
    peak = ejex(inds_dura(1)+find(integ_timeline_stim_cut(inds_dura)==max(integ_timeline_stim_cut(inds_dura)),1)-1);
else
    peak = 0;
end
aver_fr =  integ_timeline_stim_cut(find(ejex>tmin_tot,1)-1:decimate_factor:find(ejex>tmax_tot,1));
ejex=ejex(find(ejex>tmin_tot,1)-1:decimate_factor:find(ejex>tmax_tot,1));


