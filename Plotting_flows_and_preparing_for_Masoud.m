clc;
clear all;
% close all;
Coherences=[0.22 0.3 0.45 0.55];
past_window=50;
% past_time=[30 50 70 90 110 130 150 170 200]
% past_time=90;
coherence=1;
stim_resp=1;
levels=0;
figure;

past_time=30;
for coherence=1:4
if stim_resp==1
    load(['st_al_pCor_IMG_occip_front_and_Flow_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
    spans=-100:10:600;
else
    load(['rp_al_pCor_IMG_occip_front_and_Flow_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
    spans=-600:10:100;
end
%% Significance and renaming of variables

% if levels==0
%     plot(spans,nanmean(ParCorrelations_Fam_Unfam_frnt),'b')
%     hold on;
%     plot(spans,nanmean(ParCorrelations_FF_Fam_Unfam),'.-b')
%     plot(spans,nanmean(ParCorrelations_FF_Fam_Unfam_Self_Out),'-.b')
%     % %     plot(spans,nanmean(ParCorrelations_Fam_Unfam_frnt-ParCorrelations_FF_Fam_Unfam),'-c')
%     
%     plot(spans,nanmean(ParCorrelations_Fam_Unfam_ocpt),'r')
%     hold on;
%     plot(spans,nanmean(ParCorrelations_FB_Fam_Unfam),'.-r')
%     plot(spans,nanmean(ParCorrelations_FB_Fam_Unfam_Self_Out),'-.r')
%     %     plot(spans,nanmean(ParCorrelations_Fam_Unfam_ocpt-ParCorrelations_FB_Fam_Unfam),'-m')
%     %     plot(spans,[nanmean(ParCorrelations_Fam_Unfam_frnt-ParCorrelations_FF_Fam_Unfam)-nanmean(ParCorrelations_Fam_Unfam_ocpt-ParCorrelations_FB_Fam_Unfam)],'k');
% else
%     plot(spans,nanmean(ParCorrelations_Fam_Levels_frnt),'b')
%     hold on;
%     plot(spans,nanmean(ParCorrelations_FF_Fam_Levels),'.-b')
%     plot(spans,nanmean(ParCorrelations_FF_Fam_Levels_Self_Out),'-.b')
%     
%     plot(spans,nanmean(ParCorrelations_Fam_Levels_ocpt),'r')
%     hold on;
%     plot(spans,nanmean(ParCorrelations_FB_Fam_Levels),'.-r')
%     plot(spans,nanmean(ParCorrelations_FB_Fam_Levels_Self_Out),'-.r')
% end
% 

subjs=1:16;
for t=1:71
    % frontal
    
    if sum(nanmean(ParCorrelations_Fam_Unfam_frnt(:,t))>nanmean(ParCorrelations_Fam_Unfam_random_frnt(:,t,:)))>950
        signif_frnt_FamUnfam(t)=1;
    else
        signif_frnt_FamUnfam(t)=nan;
    end
    ParCorr_frnt_FamUnfam=ParCorrelations_Fam_Unfam_frnt;
    
    if sum(nanmean(ParCorrelations_FF_Fam_Unfam(:,t))>nanmean(ParCorrelations_FF_Fam_Unfam_random(:,t,:)))>950
        signif_frnt_minus_ocpt_FamUnfam(t)=1;
    else
        signif_frnt_minus_ocpt_FamUnfam(t)=nan;
    end
    ParCorr_frnt_minus_ocpt_FamUnfam=ParCorrelations_FF_Fam_Unfam;
    
    if sum(nanmean(ParCorrelations_FF_Fam_Unfam_Self_Out(:,t))>nanmean(ParCorrelations_FF_Fam_Unfam_Self_Out_random(:,t,:)))>950
        signif_frnt_minus_ocpt_and_frnt_FamUnfam(t)=1;
    else
        signif_frnt_minus_ocpt_and_frnt_FamUnfam(t)=nan;
    end
    ParCorr_frnt_minus_ocpt_and_frnt_FamUnfam=ParCorrelations_FF_Fam_Unfam_Self_Out;
    
    if sum(nanmean(ParCorrelations_Fam_Levels_frnt(:,t))>nanmean(ParCorrelations_Fam_Levels_random_frnt(:,t,:)))>950
        signif_frnt_Levels(t)=1;
    else
        signif_frnt_Levels(t)=nan;
    end
    ParCorr_frnt_Levels=ParCorrelations_Fam_Levels_frnt;
    
    if sum(nanmean(ParCorrelations_FF_Fam_Levels(:,t))>nanmean(ParCorrelations_FF_Fam_Levels_random(:,t,:)))>950
        signif_frnt_minus_ocpt_Levels(t)=1;
    else
        signif_frnt_minus_ocpt_Levels(t)=nan;
    end
    ParCorr_frnt_minus_ocpt_Levels=ParCorrelations_FF_Fam_Levels;
    
    if sum(nanmean(ParCorrelations_FF_Fam_Levels_Self_Out(:,t))>nanmean(ParCorrelations_FF_Fam_Levels_Self_Out_random(:,t,:)))>950
        signif_frnt_minus_ocpt_and_frnt_Levels(t)=1;
    else
        signif_frnt_minus_ocpt_and_frnt_Levels(t)=nan;
    end
    ParCorr_frnt_minus_ocpt_and_frnt_Levels=ParCorrelations_FF_Fam_Levels_Self_Out;
    
    % occipital
    
    if sum(nanmean(ParCorrelations_Fam_Unfam_ocpt(:,t))>nanmean(ParCorrelations_Fam_Unfam_random_ocpt(:,t,:)))>950
        signif_ocpt_FamUnfam(t)=1;
    else
        signif_ocpt_FamUnfam(t)=nan;
    end
    ParCorr_ocpt_FamUnfam=ParCorrelations_Fam_Unfam_ocpt;
    
    if sum(nanmean(ParCorrelations_FB_Fam_Unfam(:,t))>nanmean(ParCorrelations_FB_Fam_Unfam_random(:,t,:)))>950
        signif_ocpt_minus_frnt_FamUnfam(t)=1;
    else
        signif_ocpt_minus_frnt_FamUnfam(t)=nan;
    end
    ParCorr_ocpt_minus_frnt_FamUnfam=ParCorrelations_FB_Fam_Unfam;
    
    if sum(nanmean(ParCorrelations_FB_Fam_Unfam_Self_Out(:,t))>nanmean(ParCorrelations_FB_Fam_Unfam_Self_Out_random(:,t,:)))>950
        signif_ocpt_minus_frnt_and_ocpt_FamUnfam(t)=1;
    else
        signif_ocpt_minus_frnt_and_ocpt_FamUnfam(t)=nan;
    end
    ParCorr_ocpt_minus_frnt_and_ocpt_FamUnfam=ParCorrelations_FB_Fam_Unfam_Self_Out;
    
    if sum(nanmean(ParCorrelations_Fam_Levels_ocpt(:,t))>nanmean(ParCorrelations_Fam_Levels_random_ocpt(:,t,:)))>950
        signif_ocpt_Levels(t)=1;
    else
        signif_ocpt_Levels(t)=nan;
    end
    ParCorr_ocpt_Levels=ParCorrelations_Fam_Levels_ocpt;
    
    if sum(nanmean(ParCorrelations_FB_Fam_Levels(:,t))>nanmean(ParCorrelations_FB_Fam_Levels_random(:,t,:)))>950
        signif_ocpt_minus_frnt_Levels(t)=1;
    else
        signif_ocpt_minus_frnt_Levels(t)=nan;
    end
    ParCorr_ocpt_minus_frnt_Levels=ParCorrelations_FB_Fam_Levels;
    
    if sum(nanmean(ParCorrelations_FB_Fam_Levels_Self_Out(:,t))>nanmean(ParCorrelations_FB_Fam_Levels_Self_Out_random(:,t,:)))>950
        signif_ocpt_minus_frnt_and_ocpt_Levels(t)=1;
    else
        signif_ocpt_minus_frnt_and_ocpt_Levels(t)=nan;
    end
    ParCorr_ocpt_minus_frnt_and_ocpt_Levels=ParCorrelations_FB_Fam_Levels_Self_Out;
    
end




% if levels==0
%     plot(spans,signif_frnt_FamUnfam-1,'*b')
%     hold on;
%     plot(spans,signif_frnt_minus_ocpt_FamUnfam-0.99,'ob')
%     plot(spans,signif_frnt_minus_ocpt_and_frnt_FamUnfam-0.98,'.b')
%     
%     plot(spans,signif_ocpt_FamUnfam-1.01,'*r')
%     hold on;
%     plot(spans,signif_ocpt_minus_frnt_FamUnfam-1.02,'or')
%     plot(spans,signif_ocpt_minus_frnt_and_ocpt_FamUnfam-1.03,'.r')
% else
%     plot(spans,signif_frnt_Levels-1,'*b')
%     hold on;
%     plot(spans,signif_frnt_minus_ocpt_Levels-0.99,'ob')
%     plot(spans,signif_frnt_minus_ocpt_and_frnt_Levels-0.98,'.b')
%     
%     plot(spans,signif_ocpt_Levels-1.01,'*r')
%     hold on;
%     plot(spans,signif_ocpt_minus_frnt_Levels-1.02,'or')
%     plot(spans,signif_ocpt_minus_frnt_and_ocpt_Levels-1.03,'.r')
% end

for t=1:71
    if signrank(ParCorr_frnt_FamUnfam(:,t),ParCorr_frnt_minus_ocpt_FamUnfam(:,t))<0.05
        signif_frnt_1_3_FamUnfam(1,t)=1;
    else
        signif_frnt_1_3_FamUnfam(1,t)=nan;
    end
    if signrank(ParCorr_frnt_FamUnfam(:,t),ParCorr_frnt_minus_ocpt_and_frnt_FamUnfam(:,t))<0.05
        signif_frnt_1_3_FamUnfam(2,t)=1;
    else
        signif_frnt_1_3_FamUnfam(2,t)=nan;
    end
    if signrank(ParCorr_frnt_minus_ocpt_FamUnfam(:,t),ParCorr_frnt_minus_ocpt_and_frnt_FamUnfam(:,t))<0.05
        signif_frnt_1_3_FamUnfam(3,t)=1;
    else
        signif_frnt_1_3_FamUnfam(3,t)=nan;
    end
    
    
    if signrank(ParCorr_frnt_Levels(:,t),ParCorr_frnt_minus_ocpt_Levels(:,t))<0.05
        signif_frnt_1_3_Levels(1,t)=1;
    else
        signif_frnt_1_3_Levels(1,t)=nan;
    end
    if signrank(ParCorr_frnt_Levels(:,t),ParCorr_frnt_minus_ocpt_and_frnt_Levels(:,t))<0.05
        signif_frnt_1_3_Levels(2,t)=1;
    else
        signif_frnt_1_3_Levels(2,t)=nan;
    end
    if signrank(ParCorr_frnt_minus_ocpt_Levels(:,t),ParCorr_frnt_minus_ocpt_and_frnt_Levels(:,t))<0.05
        signif_frnt_1_3_Levels(3,t)=1;
    else
        signif_frnt_1_3_Levels(3,t)=nan;
    end
    
    
    if signrank(ParCorr_ocpt_FamUnfam(:,t),ParCorr_ocpt_minus_frnt_FamUnfam(:,t))<0.05
        signif_ocpt_1_3_FamUnfam(1,t)=1;
    else
        signif_ocpt_1_3_FamUnfam(1,t)=nan;
    end
    if signrank(ParCorr_ocpt_FamUnfam(:,t),ParCorr_ocpt_minus_frnt_and_ocpt_FamUnfam(:,t))<0.05
        signif_ocpt_1_3_FamUnfam(2,t)=1;
    else
        signif_ocpt_1_3_FamUnfam(2,t)=nan;
    end
    if signrank(ParCorr_ocpt_minus_frnt_FamUnfam(:,t),ParCorr_ocpt_minus_frnt_and_ocpt_FamUnfam(:,t))<0.05
        signif_ocpt_1_3_FamUnfam(3,t)=1;
    else
        signif_ocpt_1_3_FamUnfam(3,t)=nan;
    end
    
    if signrank(ParCorr_ocpt_Levels(:,t),ParCorr_ocpt_minus_frnt_Levels(:,t))<0.05
        signif_ocpt_1_3_Levels(1,t)=1;
    else
        signif_ocpt_1_3_Levels(1,t)=nan;
    end
    if signrank(ParCorr_ocpt_FamUnfam(:,t),ParCorr_ocpt_minus_frnt_and_ocpt_Levels(:,t))<0.05
        signif_ocpt_1_3_Levels(2,t)=1;
    else
        signif_ocpt_1_3_Levels(2,t)=nan;
    end
    if signrank(ParCorr_ocpt_minus_frnt_Levels(:,t),ParCorr_ocpt_minus_frnt_and_ocpt_Levels(:,t))<0.05
        signif_ocpt_1_3_Levels(3,t)=1;
    else
        signif_ocpt_1_3_Levels(3,t)=nan;
    end
end

% if levels==0
%     plot(spans,signif_frnt_1_3_FamUnfam(1,:)-1.05,'*g')
%     hold on;
%     plot(spans,signif_frnt_1_3_FamUnfam(2,:)-1.06,'*c')
%     plot(spans,signif_frnt_1_3_FamUnfam(3,:)-1.07,'*m')
%     
%     plot(spans,signif_ocpt_1_3_FamUnfam(1,:)-1.02,'og')
%     plot(spans,signif_ocpt_1_3_FamUnfam(2,:)-1.03,'oc')
%     plot(spans,signif_ocpt_1_3_FamUnfam(3,:)-1.04,'om')
% else
%     plot(spans,signif_frnt_1_3_Levels(1,:)-1.05,'*g')
%     hold on;
%     plot(spans,signif_frnt_1_3_Levels(2,:)-1.06,'*c')
%     plot(spans,signif_frnt_1_3_Levels(3,:)-1.07,'*m')
%     
%     plot(spans,signif_ocpt_1_3_Levels(1,:)-1.02,'*g')
%     plot(spans,signif_ocpt_1_3_Levels(2,:)-1.03,'*c')
%     plot(spans,signif_ocpt_1_3_Levels(3,:)-1.04,'*m')
% end

if levels==0
    difference_frontal=ParCorr_frnt_FamUnfam-ParCorr_frnt_minus_ocpt_FamUnfam;
    difference_occipital=ParCorr_ocpt_FamUnfam-ParCorr_ocpt_minus_frnt_FamUnfam;
    difference_flow_to_Frnt_minus_to_Ocpt=difference_frontal-difference_occipital;
else
    difference_frontal=ParCorr_frnt_Levels-ParCorr_frnt_minus_ocpt_Levels;
    difference_occipital=ParCorr_ocpt_Levels-ParCorr_ocpt_minus_frnt_Levels;
    difference_flow_to_Frnt_minus_to_Ocpt=difference_frontal-difference_occipital;
end
% figure;
% plot(nanmean(difference_frontal))
% hold on;
% plot(nanmean(difference_occipital))
plot(nanmean(difference_flow_to_Frnt_minus_to_Ocpt))
hold on;
end
vvv
for time=1:71
    if signrank(difference_frontal(:,time),difference_occipital(:,time))<0.05
        signif_difference_flow_to_Frnt_minus_to_Ocpt(time)=1;
    else
        signif_difference_flow_to_Frnt_minus_to_Ocpt(time)=nan;
    end
end
plot(signif_difference_flow_to_Frnt_minus_to_Ocpt*0,'*')

clearvars -except signif_difference_flow_to_Frnt_minus_to_Ocpt difference_flow_to_Frnt_minus_to_Ocpt Coherences coherence stim_resp signif_frnt_1_3_Levels signif_frnt_1_3_FamUnfam signif_ocpt_1_3_Levels signif_ocpt_1_3_FamUnfam ParCorr_frnt_minus_ocpt_and_frnt_FamUnfam ParCorr_frnt_Levels ParCorr_ocpt_FamUnfam ParCorr_ocpt_minus_frnt_and_ocpt_Levels ParCorr_ocpt_minus_frnt_Levels ParCorr_ocpt_Levels ParCorr_ocpt_minus_frnt_and_ocpt_FamUnfam ParCorr_ocpt_minus_frnt_FamUnfam ParCorr_ocpt_FamUnfam ParCorr_frnt_minus_ocpt_and_frnt_Levels ParCorr_frnt_minus_ocpt_Levels ParCorr_frnt_FamUnfam signif_frnt_Levels ParCorr_frnt_minus_ocpt_and_frnt_FamUnfam ParCorr_frnt_minus_ocpt_FamUnfam ParCorr_frnt_FamUnfam signif_frnt_FamUnfam signif_frnt_Levels signif_frnt_minus_ocpt_and_frnt_FamUnfam signif_frnt_minus_ocpt_and_frnt_Levels signif_frnt_minus_ocpt_FamUnfam signif_frnt_minus_ocpt_Levels signif_ocpt_FamUnfam signif_ocpt_Levels signif_ocpt_minus_frnt_and_ocpt_FamUnfam signif_ocpt_minus_frnt_and_ocpt_Levels signif_ocpt_minus_frnt_FamUnfam signif_ocpt_minus_frnt_Levels
% if stim_resp==1
%     save(['st_information_flow_analysis_coherence_',num2str(Coherences(coherence)),'.mat']);
% else
%     save(['rp_information_flow_analysis_coherence_',num2str(Coherences(coherence)),'.mat']);
% end
% 

