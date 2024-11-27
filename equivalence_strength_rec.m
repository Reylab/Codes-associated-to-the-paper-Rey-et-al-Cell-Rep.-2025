function equivalence_strength_rec(resp,alpha)
% resp structure as inpout and alpha = beta (type I and II errors)

twin_on=-1500;
twin_off=500;
resp_rec = 0;
            
for iresp=1:size(resp,2)
    if ~isempty(resp(iresp).trecall_phasic_ms)
        resp_rec = resp_rec + 1;
        active_cluster = resp(iresp).spike_times_Rec;
        phasic_recall = resp(iresp).trecall_phasic_ms;
        ntrials(1) = size(active_cluster{resp(iresp).responsive_Storiesindex(1)},1);
        ntrials(2) = size(active_cluster{resp(iresp).responsive_Storiesindex(2)},1);
        strength_pair = NaN*ones(max(ntrials),2);
    
        equiv_strength(resp_rec).chan = resp(iresp).channel_number;
        equiv_strength(resp_rec).class = resp(iresp).cluster;
        equiv_strength(resp_rec).pair = resp(iresp).responsive_Storiesindex;
            
        for istim=1:2 
            phas_vec = phasic_recall{resp(iresp).responsive_Storiesindex(istim)};
            spikes1 = active_cluster{resp(iresp).responsive_Storiesindex(istim)};
            spikes1 = arrayfun(@(k) spikes1{k}-phas_vec(k),[1:length(phas_vec)]','UniformOutput',false);
            strength_pair(1:ntrials(istim),istim) = cell2mat(cellfun(@(x) sum((x< twin_off) & (x> twin_on)),spikes1,'UniformOutput',0));
        end
            
    
        nsamp = min(ntrials);
        equiv_strength(resp_rec).samples = ntrials;
        equiv_strength(resp_rec).meandiff_Hz = (diff(nanmean(strength_pair)))/((twin_off-twin_on)/1000);
        equiv_strength(resp_rec).delta = sqrt(2/nsamp*(norminv(alpha)+norminv(alpha/2))^2);
        [equiv_strength(resp_rec).test_resu, equiv_strength(resp_rec).pval, equiv_strength(resp_rec).pooledSD, equiv_strength(resp_rec).meandiff] = TOST_2023(strength_pair(~isnan(strength_pair(:,1)),1), strength_pair(~isnan(strength_pair(:,2)),2), 'welch',equiv_strength(resp_rec).delta,alpha);
    end
end

pvals=cell2mat({equiv_strength.pval}');

prop_signif = sum(any(pvals>alpha,2))/size(equiv_strength,2);
fprintf('twin_on %d twin_off %d:   signif %2.1f\n',twin_on,twin_off,prop_signif*100)
