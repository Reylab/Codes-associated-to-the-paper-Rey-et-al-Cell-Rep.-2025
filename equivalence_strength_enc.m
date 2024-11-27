function equivalence_strength_enc(resp,alpha)
% resp structure as inpout and alpha = beta (type I and II errors)

twin=500;
            
for iresp=1:size(resp,2)
    active_cluster = resp(iresp).spike_times_Enc;
    ntrials = size(active_cluster{resp(iresp).responsive_Storiesindex(1)},1);
    strength_pair = zeros(ntrials,2);
    
    equiv_strength(iresp).chan = resp(iresp).channel_number;
    equiv_strength(iresp).class = resp(iresp).cluster;
    equiv_strength(iresp).pair = resp(iresp).responsive_Storiesindex;
        
    for istim=1:2 
        strength_pair(:,istim) = cell2mat(cellfun(@(x) sum((x< resp(iresp).onset_Enc(istim) + twin) & (x> resp(iresp).onset_Enc(istim))),active_cluster{resp(iresp).responsive_Storiesindex(istim)},'UniformOutput',0));
    end
        
    equiv_strength(iresp).str_diff_0 = abs(diff(median(strength_pair)));
    equiv_strength(iresp).meandiff_Hz = (diff(mean(strength_pair)))/(twin/1000);

    nsamp = min([length(strength_pair(:,1)), length(strength_pair(:,2))]);
         
    equiv_strength(iresp).delta = sqrt(2/nsamp*(norminv(alpha)+norminv(alpha/2))^2);

    [equiv_strength(iresp).test_resu, equiv_strength(iresp).pval, equiv_strength(iresp).pooledSD, equiv_strength(iresp).meandiff] = TOST_2023(strength_pair(:,1), strength_pair(:,2), 'welch',equiv_strength(iresp).delta,alpha);
end

pvals=cell2mat({equiv_strength.pval}');

prop_signif = sum(any(pvals>alpha,2))/size(equiv_strength,2);
fprintf('twin %d:   signif %2.1f\n',twin,prop_signif*100)
