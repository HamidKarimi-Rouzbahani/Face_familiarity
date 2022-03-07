clc;
clear all;
% close all;
Coherences=[0.22 0.3 0.45 0.55];
past_window=50;
% past_time=[30 50 70 90 110 130 150 170 200]
% past_time=90;
stim_resp=2;
model=1; % familiar, unfamiliar, famous, personally, self, levels
titles={'unfamiliar','familiar','levels','famous','personally','self'};
for model=1:6
    subplot(3,2,model)
    past_time=30;
    for coherence=1:4
        colors=[0.8 0.8 0.8;0.6 0.6 0.6;0.4 0.4 0.4;0.1 0.1 0.1];
        if stim_resp==1
            %             load(['st_flow_all_models_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
            load(['st_al_pCor_IMG_occip_front_and_Flow_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
            spans=-100:10:600;
        else
            %             load(['rp_flow_all_models_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
            load(['st_al_pCor_IMG_occip_front_and_Flow_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
            spans=-600:10:100;
        end
        %% Significance and renaming of variables
        if model==1
            plot(spans,nanmean(ParCorrelations_Unfamiliar_frnt-ParCorrelations_FF_Unfamiliar)-nanmean(ParCorrelations_Unfamiliar_ocpt-ParCorrelations_FB_Unfamiliar),'color',colors(coherence,:),'linewidth',3)
        elseif model==2
            plot(spans,nanmean(ParCorrelations_Familiar_frnt-ParCorrelations_FF_Familiar)-nanmean(ParCorrelations_Familiar_ocpt-ParCorrelations_FB_Familiar),'color',colors(coherence,:),'linewidth',3)
        elseif model==3
            plot(spans,nanmean(ParCorrelations_Lev_familiar_frnt-ParCorrelations_FF_Lev_familiar)-nanmean(ParCorrelations_Lev_familiar_ocpt-ParCorrelations_FB_Lev_familiar),'color',colors(coherence,:),'linewidth',3)
        elseif model==4
            plot(spans,nanmean(ParCorrelations_Famous_frnt-ParCorrelations_FF_Famous)-nanmean(ParCorrelations_Famous_ocpt-ParCorrelations_FB_Famous),'color',colors(coherence,:),'linewidth',3)
        elseif model==5
            plot(spans,nanmean(ParCorrelations_Personally_frnt-ParCorrelations_FF_Personally)-nanmean(ParCorrelations_Personally_ocpt-ParCorrelations_FB_Personally),'color',colors(coherence,:),'linewidth',3)
        elseif model==6
            plot(spans,nanmean(ParCorrelations_Self_frnt-ParCorrelations_FF_Self)-nanmean(ParCorrelations_Self_ocpt-ParCorrelations_FB_Self),'color',colors(coherence,:),'linewidth',3)
        end
        hold on;
    end
    if model==1
        legend coh=0.22 coh=0.35 coh=0.45 coh=0.55
        xlabel('time (ms)')
        ylabel('feedforward/feedback flow')
    end
    title(titles{model})
    line([spans(1) spans(end)],[0 0],'linestyle','-.','color','k')
    line([0 0],[-0.02 0.02],'linestyle','-.','color','k')
    set(gca,'FontSize',16)
end
ccc
%% old ones per coherence 

clc;
clear all;
% close all;
Coherences=[0.22 0.3 0.45 0.55];
past_window=50;
% past_time=[30 50 70 90 110 130 150 170 200]
% past_time=90;
stim_resp=2;
titles={'familiar-unfamiliar','familiarity levels'};
for model=1:2
    subplot(1,2,model)
    past_time=30;
    for coherence=1:4
        colors=[0.8 0.8 0.8;0.6 0.6 0.6;0.4 0.4 0.4;0.1 0.1 0.1];
        if stim_resp==1
            load(['st_al_pCor_IMG_occip_front_and_Flow_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
            spans=-100:10:600;
        else
            load(['rp_al_pCor_IMG_occip_front_and_Flow_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat']);
            spans=-600:10:100;
        end
        if model==1
            plot(spans,[nanmean(ParCorrelations_Fam_Unfam_frnt-ParCorrelations_FF_Fam_Unfam)-nanmean(ParCorrelations_Fam_Unfam_ocpt-ParCorrelations_FB_Fam_Unfam)],'color',colors(coherence,:),'linewidth',3);
        elseif model==2
            plot(spans,[nanmean(ParCorrelations_Fam_Levels_frnt-ParCorrelations_FF_Fam_Levels)-nanmean(ParCorrelations_Fam_Levels_ocpt-ParCorrelations_FB_Fam_Levels)],'color',colors(coherence,:),'linewidth',3);
        end
        hold on;
    end
    if model==1
        legend coh=0.22 coh=0.35 coh=0.45 coh=0.55
        xlabel('time (ms)')
        ylabel('feedforward/feedback flow')
    end
    title(titles{model})
    line([spans(1) spans(end)],[0 0],'linestyle','-.','color','k')
    line([0 0],[-0.02 0.02],'linestyle','-.','color','k')
    set(gca,'FontSize',16)
end
%% coherences pooled new ones
clc;
clear all;
% close all;
Coherences=[0.22 0.3 0.45 0.55];
past_window=50;
% past_time=[30 50 70 90 110 130 150 170 200]
% past_time=90;
stim_resp=1;
model=1; % familiar, unfamiliar, famous, personally, self, levels
titles={'unfamiliar','familiar','levels','famous','personally','self'};
for model=1:6
%     subplot(3,2,model)
    past_time=30;
    colors=ones(6,3);
    
    if stim_resp==1               
        load(['st_flow_all_models_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'.mat']);
        spans=-100:10:600;
    else
        load(['rp_flow_all_models_Novel_SP_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'.mat']);
        spans=-600:10:100;
    end
    % Significance and renaming of variables
    if model==1
        plots(model)=plot(spans,nanmean(ParCorrelations_Unfamiliar_frnt-ParCorrelations_FF_Unfamiliar)-nanmean(ParCorrelations_Unfamiliar_ocpt-ParCorrelations_FB_Unfamiliar),'linewidth',3)
    elseif model==2
        plots(model)=plot(spans,nanmean(ParCorrelations_Familiar_frnt-ParCorrelations_FF_Familiar)-nanmean(ParCorrelations_Familiar_ocpt-ParCorrelations_FB_Familiar),'linewidth',3)
    elseif model==3
        plots(model)=plot(spans,nanmean(ParCorrelations_Lev_familiar_frnt-ParCorrelations_FF_Lev_familiar)-nanmean(ParCorrelations_Lev_familiar_ocpt-ParCorrelations_FB_Lev_familiar),'linewidth',3)
    elseif model==4
        plots(model)=plot(spans,nanmean(ParCorrelations_Famous_frnt-ParCorrelations_FF_Famous)-nanmean(ParCorrelations_Famous_ocpt-ParCorrelations_FB_Famous),'linewidth',3)
    elseif model==5
        plots(model)=plot(spans,nanmean(ParCorrelations_Personally_frnt-ParCorrelations_FF_Personally)-nanmean(ParCorrelations_Personally_ocpt-ParCorrelations_FB_Personally),'linewidth',3)
    elseif model==6
        plots(model)=plot(spans,nanmean(ParCorrelations_Self_frnt-ParCorrelations_FF_Self)-nanmean(ParCorrelations_Self_ocpt-ParCorrelations_FB_Self),'linewidth',3)
    end
    hold on;
    if model==6
        legend([plots(1),plots(2),plots(3),plots(4),plots(5),plots(6)],{'unfamiliar','familiar','levels','famous','personally','self'});
        xlabel('time (ms)')
        ylabel('feedforward/feedback flow')
    end
%     title(titles{model})
    line([spans(1) spans(end)],[0 0],'linestyle','-.','color','k')
    line([0 0],[-0.01 0.01],'linestyle','-.','color','k')
    set(gca,'FontSize',16)
end
