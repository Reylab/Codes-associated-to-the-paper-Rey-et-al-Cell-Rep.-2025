function plot_all_responses_matrix(resp,phase)
% phase = 'encoding';
% phase = 'scr';
% phase = 'recall';

tmin_base = -600;
tmax_base = -100;
% smooth_bin=1500;
smooth_bin=3000;
sigma_gauss = 10;
alpha_gauss = 3.035;
half_ancho_gauss = alpha_gauss *sigma_gauss;
sample_period = 1000/30000; % sample period for the spike list - window convolution in ms/sample
decimate_factor = 1; % decimation factor used after windowing the spikes times.
N_gauss = 2*round(half_ancho_gauss/sample_period)+1; % Number of points of the gaussian window
int_window = gausswin(N_gauss, alpha_gauss);
int_window = 1000*int_window/sum(int_window)/sample_period;

% downs=30;
downs=60;

if strcmp(phase,'recall')
    tmin_epoch = -4000;
    tmax_epoch = 3000;   
else
    tmin_epoch = -900;
    tmax_epoch = 1800;
end
 
nresp = length(resp);
ejex = tmin_epoch-half_ancho_gauss:sample_period:tmax_epoch+half_ancho_gauss;  %should be the same length as spike_timeline
which_times = find(ejex>tmin_epoch,1)-1:decimate_factor:find(ejex>tmax_epoch,1);

if strcmp(phase,'scr')
    FR_mat_resp = NaN(nresp,length(which_times(1:downs:end)));
    FR_mat_non_resp = NaN(nresp,length(which_times(1:downs:end)));
    stims = cellfun(@(x, y) [x, y], {resp.responsive_stilmulus_VPindex}', {resp.nonresponsive_stilmulus_VPindex}', 'UniformOutput', false);
else
    FR_mat_resp1 = NaN(nresp,length(which_times(1:downs:end)));
    FR_mat_non_resp1 = NaN(nresp,length(which_times(1:downs:end)));
    FR_mat_resp2 = NaN(nresp,length(which_times(1:downs:end)));
    FR_mat_non_resp2 = NaN(nresp,length(which_times(1:downs:end)));
    stims = {resp.responsive_Storiesindex}';
end

latencies=cell2mat({resp.onset_Enc}');
meani=mean(latencies,2);
[~,inds_resps] = sort(meani);
resp_sorted= resp(inds_resps);

mu_base_med = NaN(nresp,1);
std_base_med= NaN(nresp,1);

% figure(1)

for iresp=1:nresp
    % clf(1)  
        
%     t_pre_ms = G.t_pre/1000;
%     t_pos_ms = G.t_pos/1000;

    if strcmp(phase,'encoding')
        active_cluster = resp(inds_resps(iresp)).spike_times_Enc;           
    elseif strcmp(phase,'recall')
        if ~isempty(resp(inds_resps(iresp)).trecall_phasic_ms)
            active_cluster = resp(inds_resps(iresp)).spike_times_Rec;
            phasic_recall = resp(inds_resps(iresp)).trecall_phasic_ms;
        else
            continue
        end
    elseif strcmp(phase,'scr')
        active_cluster = resp(inds_resps(iresp)).spike_times_VP;  
    end

    spikes{1} = active_cluster{stims{inds_resps(iresp)}(1)};
    spikes{2} = active_cluster{stims{inds_resps(iresp)}(2)};
   
    if ~strcmp(phase,'scr')
        non_resp_stims = setdiff([1:4],stims{inds_resps(iresp)});
        spikes{3} = active_cluster{non_resp_stims(1)};
        spikes{4} = active_cluster{non_resp_stims(2)};
    end

        lstim = length(active_cluster);
        ntrials_stims = zeros(lstim,1);
        sp_count_base_FR = NaN(lstim,60);
 
        for st=1:lstim
            spikes1 = active_cluster{st};
            ntrials_stims(st) = size(spikes1,1);
            if iscell(spikes1)
                sp_count_base_FR(st,1:ntrials_stims(st)) = cell2mat(cellfun(@(x) sum((x< -tmin_base) & (x> -tmax_base)),spikes1,'UniformOutput',0))';
            else
                sp_count_base_FR(st,1:ntrials_stims(st)) = sum((spikes1< -tmin_base) & (spikes1> -tmax_base),2)';
            end
        end    

        medians_base = median(sp_count_base_FR,2,'omitnan');
        mu_base_med(iresp) = mean(medians_base);
        std_base_med(iresp) = std(medians_base);

    if strcmp(phase,'recall')
    stim_inds = [stims{inds_resps(iresp)} non_resp_stims];
        for jjj=1:4
            phas_vec = phasic_recall{stim_inds(jjj)};
            spikes{jjj} = arrayfun(@(k) spikes{jjj}{k}-phas_vec(k),[1:length(phas_vec)]','UniformOutput',false);
        end
    end

    [ntrials,~] = cellfun(@size,spikes); 
    aver_fr={};

    for ii=1:length(spikes)
        if iscell(spikes{ii})
            all_spks=cell2mat(spikes{ii}');
        else
            all_spks=spikes{ii}(:);
        end

        spikes_tot=all_spks(all_spks < tmax_epoch+half_ancho_gauss & all_spks > tmin_epoch-half_ancho_gauss);
        spike_timeline = hist(spikes_tot,(tmin_epoch-half_ancho_gauss:sample_period:tmax_epoch+half_ancho_gauss))/ntrials(ii);
        n_spike_timeline = length(spike_timeline); %should be the same length as ejex
        integ_timeline_stim = conv(spike_timeline, int_window);
        integ_timeline_stim_cut = integ_timeline_stim(round(half_ancho_gauss/sample_period)+1:n_spike_timeline+round(half_ancho_gauss/sample_period));
        aver_fr{ii} =  smooth(integ_timeline_stim_cut(which_times),smooth_bin);

        % subplot(length(spikes),1,ii)
        % plot(ejex(which_times(1:downs:end)),aver_fr{ii}(1:downs:end))
        % maxi(ii) = max(aver_fr{ii});
    end        
    % sgtitle(sprintf('%d out of %d. %s_ch%d_%s',iresp,nresp,resp(inds_resps(iresp)).suffix_grapes,resp(inds_resps(iresp)).chan,resp(inds_resps(iresp)).class),'interpreter','none')
    % for ii=1:length(spikes)
    %     subplot(length(spikes),1,ii)
    %     set(gca,'ylim',[0 max(maxi)])
    % end

    if strcmp(phase,'scr')
        FR_mat_resp(iresp,:) = aver_fr{1}(1:downs:end);
        FR_mat_non_resp(iresp,:) = aver_fr{2}(1:downs:end);
    else
        FR_mat_resp1(iresp,:) = aver_fr{1}(1:downs:end);
        FR_mat_resp2(iresp,:) = aver_fr{2}(1:downs:end);
        FR_mat_non_resp1(iresp,:) = aver_fr{3}(1:downs:end);
        FR_mat_non_resp2(iresp,:) = aver_fr{4}(1:downs:end);
    end

end
ejex=ejex(which_times(1:downs:end));

if strcmp(phase,'scr')
    save(['FR_MAT_' phase],'FR_mat_resp','FR_mat_non_resp','mu_base_med','std_base_med','ejex','resp_sorted')
else
    save(['FR_MAT_' phase],'FR_mat_resp1','FR_mat_resp2','FR_mat_non_resp1','FR_mat_non_resp2','mu_base_med','std_base_med','ejex','resp_sorted')
end
    
%%
% clearvars
figure(10)
clf(10)
dir_base = pwd;
% phase = 'encoding';
% phase = 'scr';
% phase = 'recall';
which_norm = 'peakmuRenc+mu_base_sub';
with_gaps=0;
% with_gaps=1;
signif_only = 1;
% signif_only = 0;

tmin_plot = -100;
tmax_plot = 1000;
str_titles = {'Responsive 1';'Responsive 2';'Non-Responsive 1';'Non-Responsive 2'};
cuales = 1:length(resp_sorted);

if strcmp(phase,'scr')
    mats2plot{1} = FR_mat_resp;
    mats2plot{2} = FR_mat_non_resp;  
    str_titles = {'Responsive';'Non-Responsive'};
elseif strcmp(phase,'encoding')
    mats2plot{1} = FR_mat_resp1;
    mats2plot{2} = FR_mat_resp2;
    mats2plot{3} = FR_mat_non_resp1;
    mats2plot{4} = FR_mat_non_resp2;
elseif strcmp(phase,'recall')
    % submat = resp_sorted(cell2mat({resp_sorted.is_recall}));
    cuales = ~isnan(sum(FR_mat_resp1,2));
    if signif_only
        cuales =cell2mat({resp_sorted.is_signif_recall}');
    end
    mats2plot{1} = FR_mat_resp1(cuales,:);
    mats2plot{2} = FR_mat_resp2(cuales,:);
    mats2plot{3} = FR_mat_non_resp1(cuales,:);
    mats2plot{4} = FR_mat_non_resp2(cuales,:);
    tmin_plot = -3000;
    tmax_plot = 2000;
end

x = [tmin_plot tmax_plot];
colors='brkm';
mats2plot_peak ={};
leg_lines=[];
which_times = find(ejex>tmin_plot,1)-1:find(ejex>tmax_plot,1);

if strcmp(phase,'encoding')
        peak_norm_enc = zeros(length(mu_base_med(cuales)),1);
        for imat=1:2
            mat4norm = mats2plot{imat};
            if contains(which_norm,'mu_base_sub')
                mat4norm = mat4norm - repmat(mu_base_med(cuales),1,size(mat4norm,2));
                mat4norm(mat4norm<0)=0;
            end
            peak_norm_enc(:,imat)=max(mat4norm(:,which_times),[],2);
        end
        if contains(which_norm,'peakmuR')
            peak_norm_enc = mean(peak_norm_enc,2);
        elseif contains(which_norm,'peakR')
            peak_norm_enc = peak_norm_enc(:,1);
        end
%                 peak_norm_enc = max(peak_norm_enc,[],2);
        save(fullfile(dir_base,'peak_norm_enc.mat'),'peak_norm_enc')    
elseif strcmp(phase,'recall')
        peak_norm_rec = zeros(length(mu_base_med(cuales)),1);
        for imat=1:2
            mat4norm = mats2plot{imat};
            if contains(which_norm,'mu_base_sub')
                mat4norm = mat4norm - repmat(mu_base_med(cuales),1,size(mat4norm,2));
                mat4norm(mat4norm<0)=0;
            end
            peak_norm_rec(:,imat)=max(mat4norm(:,which_times),[],2);
        end
        if contains(which_norm,'peakmuR')
            peak_norm_rec = mean(peak_norm_rec,2);
        elseif contains(which_norm,'peakR')
            peak_norm_rec = peak_norm_rec(:,1);
        end        
elseif strcmp(phase,'scr')
    load(fullfile(dir_base,'peak_norm_enc.mat'))
end

for imat=1:size(mats2plot,2)
    if with_gaps
        y = [1 size(resp_sorted,2)];
    else        
        y = [1 size(mats2plot{imat},1)];
    end
    if contains(which_norm,'mu_base_sub')
        mats2plot{imat} = mats2plot{imat} - repmat(mu_base_med(cuales),1,size(mats2plot{imat},2));
        mats2plot{imat}(mats2plot{imat}<0)=0;
    end

    if ~strcmp(phase,'recall')
        mats2plot_peak{imat} = (mats2plot{imat}'*diag(1./peak_norm_enc))';
    else
        mats2plot_peak{imat} = (mats2plot{imat}'*diag(1./peak_norm_rec))';
    end

    figure(imat)
    clf(imat)
    if with_gaps
        matwgaps = zeros(size(resp_sorted,2),size(mats2plot_peak{imat},2));
        matwgaps(cell2mat({resp_sorted.is_recall}),:) = mats2plot_peak{imat};
        image(x,y,256*matwgaps(:,which_times))
    else
        image(x,y,256*mats2plot_peak{imat}(:,which_times))
    end

    xlim(x)
    colorbar
    
    title(sprintf('%s phase. %s.  %s normalization ',phase,str_titles{imat},which_norm),'interpreter','none')
    xlabel('Time from stimulus onset (ms)')
    ylabel('Response #')
    set(gca,'xminorgrid','on')
    hold on
   
    set(gca,'yticklabels',[])

    % print(imat,'-dpng',fullfile(dir_base,sprintf('%s_phase_%s_%s_norm_signonly_%d.png',phase,str_titles{imat},which_norm,signif_only)));
    % print(imat,'-depsc2','-tiff',fullfile(dir_base,sprintf('%s_phase_%s_%s_norm_signonly_%d.eps',phase,str_titles{imat},which_norm,signif_only)));

    GA= mean(mats2plot_peak{imat}(:,which_times));
    sem = std(mats2plot_peak{imat}(:,which_times))/sqrt(size(mats2plot_peak{imat},1)); 

    figure(10)
    hold on
    leg_lines(imat)=plot(ejex(which_times),GA,'Color',colors(imat),'linewidth',2);    
    jbfill(ejex(which_times),GA + sem,GA - sem,colors(imat),'none',1,0.2);
end

legend(leg_lines,str_titles,'location','best')
xlim(x)
box on
set(gca,'xminorgrid','on')
if contains(which_norm,'mu_base_sub')
    title(sprintf('%s phase. baseline subtracted ',phase),'interpreter','none')
else
    title(sprintf('%s phase. no normalization ',phase),'interpreter','none')
end
xlabel('Time from stimulus onset (ms)')
ylabel('GAvg IFR (Hz)','fontsize',12,'color','k')
%  print(10,'-dpng',fullfile(dir_base,sprintf('GAvg %s_phase_%s_norm_signonly_%d.png',phase,which_norm,signif_only)));
% print(10,'-depsc2','-tiff',fullfile(dir_base,sprintf('GAvg %s_phase_%s_norm_signonly_%d.eps',phase,which_norm,signif_only)));
set(gca,'xticklabels',[],'yticklabels',[])
%print(10,'-dsvg',fullfile(dir_base,sprintf('GAvg %s_phase_%s_norm_signonly_%d.svg',phase,which_norm,signif_only)));

