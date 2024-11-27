function surro_strength_recall(resp)

twin_on=-1500;
twin_off=500;

nsurro = 1000;
resp_rec = 0;

for iresp=1:size(resp,2)
    if ~isempty(resp(iresp).trecall_phasic_ms)
        resp_rec = resp_rec + 1;
        active_cluster = resp(iresp).spike_times_Rec;
        phasic_recall = resp(iresp).trecall_phasic_ms;
        ntrials(1) = size(active_cluster{resp(iresp).responsive_Storiesindex(1)},1);
        ntrials(2) = size(active_cluster{resp(iresp).responsive_Storiesindex(2)},1);
        strength_pair = NaN*ones(max(ntrials),2);
    
    
    surro_strength(resp_rec).chan = resp(iresp).channel_number;
    surro_strength(resp_rec).class = resp(iresp).cluster;
    surro_strength(resp_rec).pair = resp(iresp).responsive_Storiesindex;
    
    for istim=1:2
        phas_vec = phasic_recall{resp(iresp).responsive_Storiesindex(istim)};
        spikes1 = active_cluster{resp(iresp).responsive_Storiesindex(istim)};
        spikes1 = arrayfun(@(k) spikes1{k}-phas_vec(k),[1:length(phas_vec)]','UniformOutput',false);
        strength_pair(1:ntrials(istim),istim) = cell2mat(cellfun(@(x) sum((x< twin_off) & (x> twin_on)),spikes1,'UniformOutput',0));
    end
    
    surro_strength(resp_rec).str_diff_0 = abs(diff(nanmedian(strength_pair)));
    surro_strength(resp_rec).norm_str_diff_0 = abs(diff(nanmedian(strength_pair)))/max(nanmedian(strength_pair));
    strengths_noNaN = strength_pair(~isnan(strength_pair));
    for isurro=1:nsurro
        rand_ind = randperm(numel(strengths_noNaN));
        surro_strength(resp_rec).str_diff_surro(isurro) = abs(median(strengths_noNaN(rand_ind(1:ntrials(1))))-median(strengths_noNaN(rand_ind(ntrials(1)+1:end))));
    end
    surro_strength(resp_rec).pval = sum(surro_strength(resp_rec).str_diff_surro>surro_strength(resp_rec).str_diff_0)/nsurro;
    end
end

prop_signif = sum(abs(cell2mat({surro_strength.pval})-0.5)>0.45)/size(surro_strength,2);
rint=randi(1000);
fprintf('twin_on %d twin_off %d rint %d:   signif %2.1f\n',twin_on,twin_off,rint,prop_signif*100)
