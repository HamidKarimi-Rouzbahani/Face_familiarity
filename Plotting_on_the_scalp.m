clc;
clear all;
close all;
% load('st_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed_region_3_channels.mat');
% load('rp_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed_region_3_channels.mat');
% load('st_Decoding_coh_highs_lows_tr_cor_ts_cor_baselined_windowed_region_3_channels.mat')
% load('rp_Decoding_coh_highs_lows_tr_cor_ts_cor_baselined_windowed_region_3_channels.mat')

% channels={'Fp1';'Fpz';'Fp2';'F7';'F3';'Fz';'F4';'F8';'FC5';'FC1';'FC2';'FC6';'M1';'T7';'C3';'Cz';'C4';'T8';'M2';'CP5';'CP1';'CP2';'CP6';'P7';'P3';'Pz';'P4';'P8';'POz';'O1';'O2';'EOG';'AF7';'AF3';'AF4';'AF8';'F5';'F1';'F2';'F6';'FC3';'FCz';'FC4';'C5';'C1';'C2';'C6';'CP3';'CP4';'P5';'P1';'P2';'P6';'PO5';'PO3';'PO4';'PO6';'FT7';'FT8';'TP7';'TP8';'PO7';'PO8';'Oz'};
% % M1 and EOG
% X=squeeze(nanmean(accuracy(:,1,:,:),1));
% plot([1:50:2000]-1500,mean(X'));
% xlabel('Time [ms]')
% ylabel('Decoding')
% ccc
% X=vertcat(X,X);

% save('data_X_rp.mat','X');
% eeglab

%%
load('St_aligned_ERPs_categories_coh_0.22.mat')

X=squeeze(nanmean(ff(:,:,:,1),1)-nanmean(nanmean(ff(:,:,:,1),1),4));

save('X.mat','X');
eeglab