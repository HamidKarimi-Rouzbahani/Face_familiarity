clc;
clear all;
% close all;
Coherences=[0.22 0.3 0.45 0.55];
past_window=50;
% past_time=[30 50 70 90 110 130 150 170 200]
% past_time=90;
coherence=4;
stim_resp=1;
levels=0;
figure;

% for past_time=[30 70 100 150 200 250 350] %% my idea
for past_time=[30] %% my idea Pearson [30 40 50 100 ...]
    
    % for past_time=[30 70 110 150 170 200] %past window=200 stim-aligned only not
    % good
    % for past_time=[30 90 110 130 150 170 200] % past window=50 (for
    % stim/resp-aligned up to 350 and 400) best 90
    %     for past_time=[30 70 110 150 170 200 350] %past window=10 not good
    %     too
    %     if stim_resp==1
    %         load(['st_al_pCor_IMG_occip_front_and_Flow_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam_frnt','ParCorrelations_Fam_Unfam_ocpt','ParCorrelations_Fam_Levels_frnt','ParCorrelations_Fam_Levels_ocpt','ParCorrelations_Fam_Unfam_random_frnt','ParCorrelations_Fam_Levels_random_frnt','ParCorrelations_Fam_Unfam_random_ocpt','ParCorrelations_Fam_Levels_random_ocpt','ParCorrelations_FF','ParCorrelations_FB','ParCorrelations_FF_random','ParCorrelations_FB_random');
    %         spans=-100:10:600;
    %     else
    %         load(['rp_al_pCor_IMG_occip_front_and_Flow_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam_frnt','ParCorrelations_Fam_Unfam_ocpt','ParCorrelations_Fam_Levels_frnt','ParCorrelations_Fam_Levels_ocpt','ParCorrelations_Fam_Unfam_random_frnt','ParCorrelations_Fam_Levels_random_frnt','ParCorrelations_Fam_Unfam_random_ocpt','ParCorrelations_Fam_Levels_random_ocpt','ParCorrelations_FF','ParCorrelations_FB','ParCorrelations_FF_random','ParCorrelations_FB_random');
    %         spans=-600:10:100;
    %     end
    %     plot(spans,nanmean(ParCorrelations_FF)-nanmean(ParCorrelations_FB))
    %     hold on;
    %
    % plot(spans,nanmean(ParCorrelations_Fam_Unfam_frnt))
    % hold on;
    % plot(spans,nanmean(ParCorrelations_Fam_Unfam_ocpt))
    % plot(spans,nanmean(ParCorrelations_Fam_Levels_frnt))
    % plot(spans,nanmean(ParCorrelations_Fam_Levels_ocpt))
    
    
    
    %     if stim_resp==1
    %         load(['st_al_pCor_IMG_occip_front_and_Flow_Novel_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
    %         spans=-100:10:600;
    %     else
    %         load(['rp_al_pCor_IMG_occip_front_and_Flow_Novel_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
    %         spans=-600:10:100;
    %     end
    %
    % %     plot(spans,nanmean(ParCorrelations_Fam_Unfam_frnt))
    % %     hold on;
    % %     plot(spans,nanmean(ParCorrelations_FF_Fam_Unfam))
    %     plot(spans,nanmean(ParCorrelations_Fam_Unfam_ocpt))
    %     hold on;
    %     plot(spans,nanmean(ParCorrelations_FB_Fam_Unfam),'-.')
    
%         if stim_resp==1
%             load(['st_al_pCor_IMG_occip_front_and_Flow_Novel_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
%             spans=-100:10:600;
%         else
%             load(['rp_al_pCor_IMG_occip_front_and_Flow_Novel_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
%             spans=-600:10:100;
%         end
%     
%         if levels==0
%             plot(spans,nanmean(ParCorrelations_Fam_Unfam_frnt),'b')
%             hold on;
%             plot(spans,nanmean(ParCorrelations_FF_Fam_Unfam),'.-b')
%             plot(spans,nanmean(ParCorrelations_FF_Fam_Unfam_Self_Out),'-.b')
%     
%             plot(spans,nanmean(ParCorrelations_Fam_Unfam_ocpt),'r')
%             hold on;
%             plot(spans,nanmean(ParCorrelations_FB_Fam_Unfam),'.-r')
%             plot(spans,nanmean(ParCorrelations_FB_Fam_Unfam_Self_Out),'-.r')
%         else
%     
%             plot(spans,nanmean(ParCorrelations_Fam_Levels_frnt),'b')
%             hold on;
%             plot(spans,nanmean(ParCorrelations_FF_Fam_Levels),'.-b')
%             plot(spans,nanmean(ParCorrelations_FF_Fam_Levels_Self_Out),'-.b')
%     
%             plot(spans,nanmean(ParCorrelations_Fam_Levels_ocpt),'r')
%             hold on;
%             plot(spans,nanmean(ParCorrelations_FB_Fam_Levels),'.-r')
%             plot(spans,nanmean(ParCorrelations_FB_Fam_Levels_Self_Out),'-.r')
%         end
%     
    
    if stim_resp==1
        load(['st_al_pCor_IMG_occip_front_and_Flow_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
        spans=-100:10:600;
    else
        load(['rp_al_pCor_IMG_occip_front_and_Flow_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
        spans=-600:10:100;
    end
    
    if levels==0
        plot(spans,nanmean(ParCorrelations_Fam_Unfam_frnt),'b')
        hold on;
        plot(spans,nanmean(ParCorrelations_FF_Fam_Unfam),'.-b')
        plot(spans,nanmean(ParCorrelations_FF_Fam_Unfam_Self_Out),'-.b')
        plot(spans,nanmean(ParCorrelations_Fam_Unfam_frnt-ParCorrelations_FF_Fam_Unfam),'-c')

        plot(spans,nanmean(ParCorrelations_Fam_Unfam_ocpt),'r')
        hold on;
        plot(spans,nanmean(ParCorrelations_FB_Fam_Unfam),'.-r')
        plot(spans,nanmean(ParCorrelations_FB_Fam_Unfam_Self_Out),'-.r')
        plot(spans,nanmean(ParCorrelations_Fam_Unfam_ocpt-ParCorrelations_FB_Fam_Unfam),'-m')
        plot(spans,[nanmean(ParCorrelations_Fam_Unfam_frnt-ParCorrelations_FF_Fam_Unfam)-nanmean(ParCorrelations_Fam_Unfam_ocpt-ParCorrelations_FB_Fam_Unfam)],'k');
    else
        
        plot(spans,nanmean(ParCorrelations_Fam_Levels_frnt),'b')
        hold on;
        plot(spans,nanmean(ParCorrelations_FF_Fam_Levels),'.-b')
        plot(spans,nanmean(ParCorrelations_FF_Fam_Levels_Self_Out),'-.b')
        
        plot(spans,nanmean(ParCorrelations_Fam_Levels_ocpt),'r')
        hold on;
        plot(spans,nanmean(ParCorrelations_FB_Fam_Levels),'.-r')
        plot(spans,nanmean(ParCorrelations_FB_Fam_Levels_Self_Out),'-.r')
    end
    
    
end


