function surro_spklat_enc(resp)

par.numper = 1000;

sigma_gauss = 10;
par.nstd = 4;
par.nstd_fact = 4;
par.alpha_gauss = 3.035;
par.half_ancho_gauss = par.alpha_gauss * sigma_gauss;
par.lat_scale = 1;
par.thr_min = 5;
par.over_threshold_time = 90;
par.below_threshold_time = 100;
par.t_down=20;
par.tmin_base = -900;
par.spk_ons_min = 100; % start searching for onset
par.tmin_baseFR = par.tmin_base;
par.tmax_base = -100;
par.tmax_post = 1500;
par.nanin = 1; 
par.nanthres=0;

for iresp=1:size(resp,2)
    active_cluster = resp(iresp).spike_times_Enc;
 
    surro_spklat(iresp).chan = resp(iresp).channel_number;
    surro_spklat(iresp).class = resp(iresp).cluster;
    surro_spklat(iresp).pair = resp(iresp).responsive_Storiesindex;
   
    lstim = length(active_cluster);
    ntrials_stims = zeros(lstim,1);
    all_spks = [];
    for st=1:lstim
        spikes1 = active_cluster{st};
        ntrials_stims(st) = size(spikes1,1);
        all_spks=[all_spks cell2mat(spikes1')];
    end
    base_spikes = all_spks((all_spks< par.tmax_base+par.half_ancho_gauss) & (all_spks> par.tmin_baseFR-par.half_ancho_gauss));
    num_base = sum(ntrials_stims);
    
    nspks1 = cell2mat(cellfun(@length,active_cluster{resp(iresp).responsive_Storiesindex(1)},'UniformOutput',0));
    nspks2 = cell2mat(cellfun(@length,active_cluster{resp(iresp).responsive_Storiesindex(2)},'UniformOutput',0));
    
    spikes1 = NaN*ones(ntrials_stims(resp(iresp).responsive_Storiesindex(1)),max([max(nspks1) max(nspks2)]));
    for jj=1:ntrials_stims(resp(iresp).responsive_Storiesindex(1))
        spikes1(jj,1:nspks1(jj)) = active_cluster{resp(iresp).responsive_Storiesindex(1)}{jj};
    end    
    ind_window =(spikes1< par.tmax_post+par.half_ancho_gauss) & (spikes1> min(par.tmin_base,-1000)-par.half_ancho_gauss);
    spikes1(~ind_window)=nan;
    sp1=spikes1;
    
    spikes1 = NaN*ones(ntrials_stims(resp(iresp).responsive_Storiesindex(2)),max([max(nspks1) max(nspks2)]));
    for jj=1:ntrials_stims(resp(iresp).responsive_Storiesindex(2))
        spikes1(jj,1:nspks2(jj)) = active_cluster{resp(iresp).responsive_Storiesindex(2)}{jj};
    end    
    ind_window =(spikes1< par.tmax_post+par.half_ancho_gauss) & (spikes1> min(par.tmin_base,-1000)-par.half_ancho_gauss);
    spikes1(~ind_window)=nan;
    sp2=spikes1;
    
    fprintf('%d/%d ... ',iresp,size(resp,2))
    [p,p_ok,numnan,~,latency1,latency2,realdiff,permdiff] = latency_test(sp1,sp2,base_spikes,num_base,par);
        
    fprintf('DONE \n')

    surro_spklat(iresp).latency1 = latency1;
    surro_spklat(iresp).latency2 = latency2;
    surro_spklat(iresp).realdiff = realdiff;
    surro_spklat(iresp).permdiff = permdiff;
    surro_spklat(iresp).p = p;
    surro_spklat(iresp).p_ok = p_ok;
    surro_spklat(iresp).numnan = numnan;   
    
end

prop_signif = sum(cell2mat({surro_spklat.p})<0.05)/size(surro_spklat,2);
fprintf('twin %d rint %d:   signif %2.1f\n',twin,rint,prop_signif*100)
