function [p,p_ok,numnan,numnanok,lat1,lat2,realdiff,permdiff,p_dura,numnan_dura,dura1,dura2,realdiff_dura,permdiff_dura] = latency_test(sp1,sp2,base_spikes,num_base,par)
% Permutation test on two populations of spikes (sp1 and sp2) to test the hypothesis (H0)
%that the two populations have the same latency.

%INPUT:   sp1 (nt,xx)
%         sp2 (nt,xx)       if par.kind==gauss, these are matrices containing the spikes time in the response window organized
%                           by trial, the two populations can have different number of trials all the empty spaces in the matrices are filled by nans
%                           if par.kind==poisson, these are the single trial latencies (column 1) and offset (column2),
%                           (nan if not estimated)
%         base_spikes       all the spikes'times in baseline for all the stimuli (empty for kind==poisson)
%         num_base          the total number of trials (should be
%                           trials*number of stimuli) from which the base_spikes were
%                           obtained (empty for kind==poisson)
%         par               contains all the parameters to compute the latency



%OUTPUT:  p                 is the significance level to reject the hypothesis (i.e. if p=< 0.01 the Ho is
%                           rejected and the two latencies can be considered different at a significance level p
%         lat1 and lat2:    latency for sp1 and sp2 respectively
%         realdiff          absolute difference betwee lat1 and lat2
%         permdiff(numper,2):the first column contains the absolute
%                            latency diffrence for each permutation,
%                           the second column indicates whether the latency was valid or not:
%                              1   Both latencies were valid (isresp=='y' for both)
%                              0   at least one of the two latencies was not a resp (it passed the threshold, but the duration wasn't long enough)
%                             -1  at least one of the two latencies was not a resp (it didn't pass the threshold in any point)

if ~isfield(par,'numper') || isempty(par.numper), par.numper=999; end
if ~isfield(par,'kind') || isempty(par.kind), par.kind='gauss'; end

nt1=size(sp1,1); %number of trials
nt2=size(sp2,1);
nt=nt1+nt2;
spk=[sp1; sp2];


if strcmp(par.kind,'gauss')
    whole1=sp1(~isnan(sp1)); whole1=whole1(:);
    whole2=sp2(~isnan(sp2)); whole2=whole2(:);
    [lat1,~,dura1,~]=compute_latency(whole1,base_spikes,num_base,nt1,par);
    [lat2,~,dura2,~]=compute_latency(whole2,base_spikes,num_base,nt2,par);
elseif strcmp(par.kind,'poisson')   
    lat1 = nanmedian(sp1(:,1));
    lat2 = nanmedian(sp2(:,1));
    dura1=nanmedian(sp1(:,2))-lat1;
    dura2=nanmedian(sp2(:,2))-lat2;
end

if isnan(lat1) || isnan(lat2)
    p = nan; p_ok = nan; numnan = nan; numnanok = nan; lat1=nan; lat2=nan; realdiff=nan; permdiff = nan;
    p_dura = nan; numnan_dura =nan; dura1 = nan; dura2 = nan; realdiff_dura =nan; permdiff_dura = nan;
else      
    % real difference between latencies
    realdiff= abs(lat1-lat2);
    realdiff_dura= abs(dura1-dura2);
    
    %permutation test
    permdiff= nan(par.numper,2);
    latper1=nan(par.numper,2);
    latper2=nan(par.numper,2);
      
    permdiff_dura= nan(par.numper,2);
    duraper1=nan(par.numper,1);
    duraper2=nan(par.numper,1);
    for i=1:par.numper
        clear whole_perm1 whole_perm2
        if i==round(par.numper/3)
            fprintf('%d/%d done... ',i,par.numper);
        elseif i==round(2*par.numper/3)
            fprintf('%d/%d done... ',i,par.numper);
        end
        
        perm = randperm(nt);  %permute the trials
        group1= spk(perm(1:nt1),:);
        group2= spk(perm(nt1+1:nt),:);
        
        if strcmp(par.kind,'gauss')
            whole_perm1=group1(~isnan(group1)); whole_perm1=whole_perm1(:);
            whole_perm2=group2(~isnan(group2)); whole_perm2=whole_perm2(:);
            %compute permuted latencies
            [latper1(i,1),latper1(i,2),duraper1(i,1),~]=compute_latency(whole_perm1,base_spikes,num_base,nt1,par);
            [latper2(i,1),latper2(i,2),duraper2(i,1),~]=compute_latency(whole_perm2,base_spikes,num_base,nt2,par);
        elseif strcmp(par.kind,'poisson')           
            latper1(i,1) = nanmedian(group1(:,1));
            latper2(i,1) = nanmedian(group2(:,1));
            latper1(i,2) = 1; %always possible to estimate latency with poisson
            latper2(i,2) = 1;
            duraper1(i,1)=nanmedian(group1(:,2))-latper1(i,1);
            duraper2(i,1)=nanmedian(group2(:,2))-latper2(i,1);
        end 
        
        permdiff(i,1)= abs(latper1(i,1)-latper2(i,1));
        permdiff_dura(i,1)=abs(duraper1(i)-duraper2(i));
        if latper1(i,2)==-1 || latper2(i,2)==-1
            permdiff(i,2)=-1;
            permdiff_dura(i,2)=-1;
        elseif latper1(i,2)==0 || latper2(i,2)==0
            permdiff(i,2)=0;
            permdiff_dura(i,2)=0;
        else
            permdiff(i,2)=1;
            permdiff_dura(i,2)=1;
        end
    end
    nanind=isnan(permdiff(:,1));
    nanindok=isnan(permdiff(:,1)) | permdiff(:,2)==0;
    numnan=sum(nanind);
    numnanok=sum(nanindok);
    if strcmp(par.kind,'gauss')
        if par.nanin
            sorted = [permdiff(nanind); sort([permdiff(~nanind,1); realdiff])];
            sorted_ok = [ nan(numnanok,1); sort([permdiff(~nanindok,1); realdiff])];
        else
            sorted = [sort([permdiff(~nanind,1); realdiff])];
            sorted_ok = [sort([permdiff(~nanindok,1); realdiff])];
        end
    elseif strcmp(par.kind,'poisson')
        sorted = [sort([permdiff(~nanind,1); realdiff])];
    end
    rank=find(sorted>=realdiff);
    rank_ok=find(sorted_ok>=realdiff);
     if numnan/par.numper > par.nanthres || par.nanthres==0
        p = (length(sorted)+1-(rank(1)))/(length(sorted));
    else
        p=nan;
    end
    if numnanok/par.numper < par.nanthres || par.nanthres==0
        p_ok = (length(sorted_ok)+1-(rank_ok(1)))/(length(sorted_ok));
    else
        p_ok=nan;
    end
    %duration
    nanind_dura=isnan(permdiff_dura(:,1));
    numnan_dura=sum(nanind_dura);
    if strcmp(par.kind,'gauss')
        sorted = [permdiff_dura(nanind_dura); sort([permdiff_dura(~nanind_dura,1); realdiff_dura])];
        rank=find(sorted>=realdiff_dura);
        p_dura = (par.numper+2-(rank(1)))/(par.numper+1);
    elseif strcmp(par.kind,'poisson')
        if numnan_dura<(length(permdiff_dura(:,1))/2) && ~isnan(realdiff_dura) %delete the nans in the poisson case
            sorted = [sort([permdiff_dura(~nanind_dura,1); realdiff_dura])];
            rank=find(sorted>=realdiff_dura);
            p_dura = (length(sorted)+1-(rank(1)))/(length(sorted));
        else
            p_dura=nan;
        end
    end
end


