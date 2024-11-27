function [onset,is_resp_code,dura,is_resp] = compute_latency(whole,base_spikes,num_base,num_ok,par)
%if poisson method is chosen, whole input has to be a structure with the
%spikes by trial
is_resp_code=1;
is_resp = 'y';
if ~isfield(par,'kind') || isempty(par.kind), par.kind='gauss'; end

if strcmp(par.kind,'gauss')
    if ~isempty(whole)
        [overthreshold_interval_first,ups,durations,~,~,~,~, threshold, thr_orig, mu_base, aver_fr,ejex_new] = get_latency_newest_ons(whole,base_spikes,num_base,par.alpha_gauss,30000,par.half_ancho_gauss,par.tmin_baseFR,par.tmax_base,par.spk_ons_min,par.tmin_base,par.tmax_post,num_ok,par.nstd,par.nstd_fact,par.over_threshold_time,par.t_down,par.lat_scale,par.thr_min);
        
        if isempty(overthreshold_interval_first)
            [~,ss]=max(durations);
            if durations==-1
                onset = NaN;
                dura = NaN;
                is_resp_code=-1;
            else
                onset = ups(ss);
                dura = durations(ss);
                is_resp_code=0;
            end
            is_resp = 'n';
        else
            onset = ups(overthreshold_interval_first);
            off_ind = overthreshold_interval_first;
            cc=overthreshold_interval_first+1;
            while cc<=length(ups)
                if ups(cc) - (durations(cc-1)+ups(cc-1)) < par.below_threshold_time
                    off_ind = off_ind + 1;
                    cc = cc+1;
                else
                    cc=length(ups)+1;
                end
            end
            dura = durations(off_ind) + ups(off_ind) - onset;
        end
    else
        is_resp ='empty';
        is_resp_code=-2;
        onset=nan;
        dura=nan;
    end
elseif strcmp(par.kind,'poisson')
     [overthreshold_interval_first,~,~,~,~,~,~,~, ~, mu_base, ~,~] = get_latency_newest_ons(cell2mat(whole),base_spikes,num_base,par.alpha_gauss,30000,par.half_ancho_gauss,par.tmin_baseFR,par.tmax_base,par.spk_ons_min,par.tmin_base,par.tmax_post,num_ok,par.nstd,par.nstd_fact,par.over_threshold_time,par.t_down,par.lat_scale,par.thr_min);
     spk_ons = NaN*ones(num_ok,1);
     spk_off = NaN*ones(num_ok,1);

    for h=1:num_ok
        BOB=[]; EOB=[];
        t=whole{h};
        if ~isempty(t)
            t_resp = t(t>par.spk_ons_min & t<par.tmax_post_plot);
            if ~isempty(t_resp)
                if mu_base<par.base_min_thr || numel(t_resp)==1
                    spk_ons(h)=t_resp(1);
                else
                    [BOB, EOB, ~]=p_burst(t_resp, mu_base/1000);
                    if ~isempty(BOB)
                        spk_ons(h)=t_resp(BOB(1));
                        spk_off(h)=t_resp(EOB(1));
                    end
                end
            end
        end
    end
    spk_ons_burst = nanmedian(spk_ons);
    spk_off_burst = nanmedian(spk_off);
    spk_ons_burst_corrected = spk_ons;
    spk_ons_burst_corrected(abs(spk_ons-spk_ons_burst)>par.ons_margin)=NaN;
    spk_ons_burst_corrected = nanmedian(spk_ons_burst_corrected);
    
    onset=spk_ons_burst;
    dura=spk_off_burst-onset;
    if isempty(overthreshold_interval_first)
        is_resp='n';
    end
    if isnan(spk_ons_burst)
        is_resp_code=-1;
    end
end
end