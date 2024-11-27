function surro_strength_enc(resp)

twin=500;

nsurro = 1000;
            
for iresp=1:size(resp,2)
    active_cluster = resp(iresp).spike_times_Enc;
    ntrials = size(active_cluster{resp(iresp).responsive_Storiesindex(1)},1);
    strength_pair = zeros(ntrials,2);
    
    surro_strength(iresp).chan = resp(iresp).channel_number;
    surro_strength(iresp).class = resp(iresp).cluster;
    surro_strength(iresp).pair = resp(iresp).responsive_Storiesindex;

    for istim=1:2 
        strength_pair(:,istim) = cell2mat(cellfun(@(x) sum((x< resp(iresp).onset_Enc(istim) + twin) & (x> resp(iresp).onset_Enc(istim))),active_cluster{resp(iresp).responsive_Storiesindex(istim)},'UniformOutput',0));
    end
                        
    surro_strength(iresp).str_diff_0 = abs(diff(median(strength_pair)));
    surro_strength(iresp).norm_str_diff_0 = abs(diff(median(strength_pair)))/max(median(strength_pair));
    
    for isurro=1:nsurro
        rand_ind = randperm(numel(strength_pair));
        surro_strength(iresp).str_diff_surro(isurro) = abs(median(strength_pair(rand_ind(1:ntrials)))-median(strength_pair(rand_ind(ntrials+1:end))));
    end
    surro_strength(iresp).pval = sum(surro_strength(iresp).str_diff_surro>surro_strength(iresp).str_diff_0)/nsurro;
end

prop_signif = sum(abs(cell2mat({surro_strength.pval})-0.5)>0.45)/size(surro_strength,2);
rint=randi(1000);
fprintf('twin %d rint %d:   signif %2.1f\n',twin,rint,prop_signif*100)
