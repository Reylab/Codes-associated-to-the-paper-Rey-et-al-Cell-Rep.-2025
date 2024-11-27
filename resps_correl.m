close all
nsurro=10000;

for gg=1:7
    figure(gg)
    set(gg,'Color',[1 1 1], 'PaperUnits', 'inches', 'PaperType', 'A4', 'PaperPositionMode', 'auto','Units','normalized', 'OuterPosition',[0 0 1 1]);
end

load FR_MAT_scr;
tini_enc = 100;
tend_enc = 800;
tini_rec = -2500;
tend_rec = 750;

scrR = mean(FR_mat_resp(:,find(ejex>tini_enc,1,'first'):find(ejex>tend_enc,1,'first'))');
scrNR = mean(FR_mat_non_resp(:,find(ejex>tini_enc,1,'first'):find(ejex>tend_enc,1,'first'))');

scr = scrR - scrNR;  % resp - non_resp

load FR_MAT_encoding;

enc1R = mean(FR_mat_resp1(:,find(ejex>tini_enc,1,'first'):find(ejex>tend_enc,1,'first'))');
enc1NR = mean(FR_mat_non_resp1(:,find(ejex>tini_enc,1,'first'):find(ejex>tend_enc,1,'first'))');
enc2R = mean(FR_mat_resp2(:,find(ejex>tini_enc,1,'first'):find(ejex>tend_enc,1,'first'))');
enc2NR = mean(FR_mat_non_resp2(:,find(ejex>tini_enc,1,'first'):find(ejex>tend_enc,1,'first'))');

enc1 = enc1R - enc1NR;  % resp - non_resp
enc2 = enc2R - enc2NR;  % resp - non_resp

%within neurons vs between neurons
dist_med = mean((enc1-enc2).^2); 
dist_med_surr =zeros(1,nsurro);
for i=1:nsurro
    enc2_surr = enc2(randperm(numel(enc2)));
    dist_med_surr(i) = mean((enc1-enc2_surr).^2);
end
pval_enc = sum(dist_med_surr<dist_med)/nsurro;

[rho_enc,prho_enc]=corr(enc1',enc2');
%encoding 1 vs. encoding 2
figure(1)
plot(enc1,enc2,'*')
axis([1 10 2 11])
axis square
xlabel('enc1')
ylabel('enc2')
box on
title(sprintf('Encoding 1 vs. Encoding 2. rho = %1.2f, p = %1.1e. psurro= %1.1e',rho_enc,prho_enc,pval_enc))

%within neurons vs between neurons
mu_enc = mean([enc1;enc2]);
dist_med = mean((mu_enc-scr).^2); 
dist_med_surr =zeros(1,nsurro);
for i=1:nsurro
    scr_surr = scr(randperm(numel(scr)));
    dist_med_surr(i) = mean((mu_enc-scr_surr).^2);
end
pval_enc_scr = sum(dist_med_surr<dist_med)/nsurro;

[rho_enc_scr,prho_enc_scr]=corr(mean([enc1;enc2])',scr');
[rho_enc_scr_rep,prho_enc_scr_rep]=corr([enc1 enc2]',[scr scr]');
%encoding vs screening
figure(2)
line([scr(:) scr(:)]',[enc1(:) enc2(:)]','linewidth',2,'marker','*')
axis([0 11 0 11])
axis square
xlabel('scr')
ylabel('enc')
box on
title(sprintf('Screening vs. Encoding 1&2. rho = %1.2f, p = %1.1e. psurro= %1.1e',rho_enc_scr,prho_enc_scr,pval_enc_scr))

%RECALL
load FR_MAT_recall
% is_recall = [1 1 1 0 1 1 0 0 1 0 1 1 1 0 0 1 1 1 1 0 1 1 0 0 1 0 1 0 1 1 1 1 1];
% is_signif = [1 1 1 0 1 0 0 0 1 0 1 0 1 0 0 1 1 1 1 0 1 1 0 0 1 0 0 0 0 1 0 0 1];
% boolean_sign = ismember(find(is_recall),find(is_signif));
boolean_sign =cell2mat({resp_sorted.is_signif_recall}');

rec1R = mean(FR_mat_resp1(:,find(ejex>tini_rec,1,'first'):find(ejex>tend_rec,1,'first'))');
rec1NR = mean(FR_mat_non_resp1(:,find(ejex>tini_rec,1,'first'):find(ejex>tend_rec,1,'first'))');
rec2R = mean(FR_mat_resp2(:,find(ejex>tini_rec,1,'first'):find(ejex>tend_rec,1,'first'))');
rec2NR = mean(FR_mat_non_resp2(:,find(ejex>tini_rec,1,'first'):find(ejex>tend_rec,1,'first'))');

rec1 = rec1R - rec1NR;  % resp - non_resp
rec2 = rec2R - rec2NR;  % resp - non_resp
is_recall = ~isnan(rec1);

%within neurons vs between neurons
dist_med = mean((rec1-rec2).^2,'omitnan'); 
dist_med_surr =zeros(1,nsurro);
for i=1:nsurro
    rec2_surr = rec2(randperm(numel(rec2)));
    dist_med_surr(i) = mean((rec1-rec2_surr).^2,'omitnan');
end
pval_rec = sum(dist_med_surr<dist_med)/nsurro;

[rho_rec,prho_rec]=corr(rec1(is_recall)',rec2(is_recall)');
[rho_rec_sign,prho_rec_sign]=corr(rec1(boolean_sign)',rec2(boolean_sign)');
%recall 1 vs. recall 2
figure(3)
hold on
plot(rec1(boolean_sign),rec2(boolean_sign),'r*')
axis([0 9 0 9])
axis square
xlabel('rec1')
ylabel('rec2')
box on
title(sprintf('Recall 1 vs. Recall 2. rho = %1.2f, p = %1.1e; rho_sig = %1.2f, p_sig = %1.1e. psurro= %1.1e',rho_rec,prho_rec,rho_rec_sign,prho_rec_sign,pval_rec))

%screening vs recall
scr_rec_sign = scr(boolean_sign);

%within neurons vs between neurons
mu_rec = mean([rec1(boolean_sign);rec2(boolean_sign)]);
dist_med = mean((mu_rec-scr_rec_sign).^2); 
dist_med_surr =zeros(1,nsurro);
for i=1:nsurro
    scr_surr = scr_rec_sign(randperm(numel(scr_rec_sign)));
    dist_med_surr(i) = mean((mu_rec-scr_surr).^2);
end
pval_rec_scr = sum(dist_med_surr<dist_med)/nsurro;

[rho_scr_rec_sign,prho_scr_rec_sign]=corr(mean([rec1(boolean_sign);rec2(boolean_sign)])',scr_rec_sign');
figure(5)
line([scr_rec_sign(:) scr_rec_sign(:)]',[rec1(boolean_sign)' rec2(boolean_sign)']','linewidth',2,'marker','*')
axis([1.7 11.8 0 10.1])
axis square
xlabel('scr')
ylabel('rec')
box on
title(sprintf('Screening vs. Recall 1&2 SIGN ONLY. rho_sig = %1.2f, p_sig = %1.1e. psurro= %1.1e',rho_scr_rec_sign,prho_scr_rec_sign,pval_rec_scr))

%%
%encoding vs recall

enc1_rec_sign = enc1(boolean_sign);
enc2_rec_sign = enc2(boolean_sign);

%within neurons vs between neurons
all_rec = [rec1(boolean_sign) rec2(boolean_sign)];
all_enc = [enc1_rec_sign enc2_rec_sign];
dist_med = mean((all_rec-all_enc).^2); 
dist_med_surr =zeros(1,nsurro);
for i=1:nsurro
    enc_surr = all_enc(randperm(numel(all_enc)));
    dist_med_surr(i) = mean((all_rec-enc_surr).^2);
end
pval_rec_enc = sum(dist_med_surr<dist_med)/nsurro;

[rho_enc_rec_sign,prho_enc_rec_sign]=corr(all_rec',all_enc');
[rho_enc_rec_sign_mu,prho_enc_rec_sign_mu]=corr(mean([rec1(boolean_sign);rec2(boolean_sign)])',mean([enc1_rec_sign;enc2_rec_sign])');
figure(7)
line([enc1_rec_sign(:) enc2_rec_sign(:)]',[rec1(boolean_sign)' rec2(boolean_sign)']','linewidth',2,'marker','*')
axis([1.5 10.5 0 9])
axis square
xlabel('enc')
ylabel('rec')
box on
title(sprintf('Encoding 1&2 vs. Recall 1&2 SIGN ONLY. rho_sig = %1.2f, p_sig = %1.1e. psurro= %1.1e',rho_enc_rec_sign,prho_enc_rec_sign,pval_rec_enc))
