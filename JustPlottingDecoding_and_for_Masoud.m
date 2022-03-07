clc;
% close all;
clear all;

trials=2; % all trials =1; train_cor-test_cor =2; train_cor-test_incor =3; traineed_on_0.55-test_on_0.22 =4; 
          % one coherence level against all others =5;
levels=1; % for each level of coherence=1; pooled across coherence levels= 0;
region=3;
stim_resp=2;
CatCoh=1; %categories=1; coherences=2

coherences=[0.22 0.3 0.45 0.55];
steps=10;
if stim_resp==1
    span=[-500:steps:1490];
else
    span=[-1500:steps:490];
end
figure;
if trials==1
    if levels==1
        for coherence=1:4
            if stim_resp==1
                load(['st_Decoding_fam_unfam_tr_all_ts_all_baselined_windowed_region_',num2str(region),'_Coh',num2str(coherences(coherence)),'V2.mat']);
                accuracy_coherences_st(coherence,:,:,:)=accuracy;
                
            else
                load(['rp_Decoding_fam_unfam_tr_all_ts_all_baselined_windowed_region_',num2str(region),'_Coh',num2str(coherences(coherence)),'V2.mat']);
                accuracy_coherences_rp(coherence,:,:,:)=accuracy;
                
            end
            %             figure;
            
            for i=1:1
                plot(span,smooth(squeeze(nanmean(accuracy(:,i,:),1)),5),'linewidth',3);
                hold on;
            end
        end
        line([span(1) span(end)],[0.5 0.5])
        line([0 0],[0.4 0.65])
        grid on;
        %         legend ('Familiar vs Unfamiliar','Unfamiliar','Familiar','Famous','Personally familiar','Self','Location','southwest');
        legend ('Coherence = 0.22','Coherence = 0.30','Coherence = 0.45','Coherence = 0.55','Location','northwest');
        xlabel('Time [ms]')
        ylabel('Decoding accuracy')
    else
      if CatCoh==1
        if stim_resp==1
            load(['st_Decoding_fam_unfam_tr_all_ts_all_baselined_windowed_region_',num2str(region),'.mat']);
        else
            load(['rp_Decoding_fam_unfam_tr_all_ts_all_baselined_windowed_region_',num2str(region),'.mat']);
        end
      else          
        if stim_resp==1
            load(['st_Decoding_coh_low_high_tr_all_ts_all_baselined_windowed_region_',num2str(region),'.mat']);
        else
            load(['rp_Decoding_coh_low_high_tr_all_ts_all_baselined_windowed_region_',num2str(region),'.mat']);
        end
      end
    end
elseif trials==3
    if CatCoh==1
        if stim_resp==1
            load(['st_Decoding_fam_unfam_tr_cor_ts_incor_baselined_windowed_region_',num2str(region),'_New3.mat']);
        else
            load(['rp_Decoding_fam_unfam_tr_cor_ts_incor_baselined_windowed_region_',num2str(region),'_New3.mat']);
        end
    else
        if stim_resp==1
            load(['st_Decoding_coh_low_high_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat'],'accuracy');
        else
            load(['rp_Decoding_coh_low_high_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat'],'accuracy');
        end
    end
elseif trials==2
    if levels==0
        if stim_resp==1
            load(['st_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed_region_',num2str(region),'_v1.mat']);
        else
            load(['rp_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed_region_',num2str(region),'_v1.mat']);
        end
    else
        for coherence=1:4
            if stim_resp==1
                load(['st_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed_region_',num2str(region),'_Coh',num2str(coherences(coherence)),'V2.mat']);
            else
                load(['rp_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed_region_',num2str(region),'_Coh',num2str(coherences(coherence)),'V2.mat']);
            end
            %             figure;
            for i=1:1
                plot(span,smooth(squeeze(nanmean(accuracy(:,i,:),1)),5),'linewidth',3);
                hold on;
            end
        end
        
        line([span(1) span(end)],[0.5 0.5])
        line([0 0],[0.4 0.65])
        %         legend ('Familiar vs Unfamiliar','Unfamiliar','Familiar','Famous','Personally familiar','Self','Location','southwest');
        grid on;
        legend ('Coherence = 0.22','Coherence = 0.30','Coherence = 0.45','Coherence = 0.55','Location','northwest');
        xlabel('Time [ms]')
        ylabel('Decoding accuracy')
        error('What you see is wrong. If you need it I will run it!!!')
    end
    
elseif trials==4
    if stim_resp==1
        load(['st_Decoding_fam_unfam_tr_0.55_ts_0.22_region_',num2str(region),'.mat'],'accuracy');
    else
        load(['rp_Decoding_fam_unfam_tr_0.55_ts_0.22_region_',num2str(region),'.mat'],'accuracy');
    end
elseif trials==5
    if stim_resp==1
        load(['st_Decoding_coh_low_higher_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat']);
    else
        load(['rp_Decoding_coh_low_higher_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat']);
    end
    acc1=accuracy;
    if stim_resp==1
        load(['st_Decoding_coh_low_higher2_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat']);
    else
        load(['rp_Decoding_coh_low_higher2_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat']);
    end
    acc2=accuracy;
    if stim_resp==1
        load(['st_Decoding_coh_low_higher3_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat']);
    else
        load(['rp_Decoding_coh_low_higher3_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat']);
    end
    acc3=accuracy;
    if stim_resp==1
        load(['st_Decoding_coh_low_higher4_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat']);
    else
        load(['rp_Decoding_coh_low_higher4_tr_cor_ts_inc_baselined_windowed_region_',num2str(region),'.mat']);
    end
    acc4=accuracy;
    accuracy=(acc1+acc2+acc3+acc4)./4;
end


if (levels==0 && trials==1)|| (levels==0 && trials==2) || trials==3 || trials==4 || trials==5
    for i=1:6
        plot(span,smooth(squeeze(nanmean(accuracy(:,i,:),1)),5),'linewidth',3);
        hold on;
    end
    line([span(1) span(end)],[0.5 0.5])
    line([0 0],[0.4 0.8])
    grid on;
    if CatCoh==1
        legend ('Familiar vs Unfamiliar','Unfamiliar','Familiar','Famous','Personally familiar','Self','Location','northwest');
    else
        if trials~=5
            legend ('Coh=0.22 vs Coh=0.55','Coh=0.22','Coh=0.55','Location','southeast');
        else
            legend ('Coh decoding','Coh=main','Coh=other 3','Location','southeast');
        end
    end
    xlabel('Time [ms]')
    ylabel('Decoding accuracy')
end

accuracy_coherences_st=accuracy;
ccc
%% significacne analysis for Masoud
subjects=18;
figure;
for coherence=1:4
    for decoding_type=1:6
        for time=1:200
            if stim_resp==1
                p(coherence,decoding_type,time)=signrank(zeros(subjects,1)+0.5,accuracy_coherences_st(coherence,~isnan(accuracy_coherences_st(coherence,:,decoding_type,time)),decoding_type,time)');
            else
                p(coherence,decoding_type,time)=signrank(zeros(subjects,1)+0.5,accuracy_coherences_rp(coherence,~isnan(accuracy_coherences_rp(coherence,:,decoding_type,time)),decoding_type,time)');
            end
            
        end
    end
end
if stim_resp==1
    Corrected_p_values_st=ones(4,6,200);
    for coherence=1:4
        for decoding_type=1:6
            Corrected_p_values_st(coherence,decoding_type,40:110)=(squeeze(p(coherence,decoding_type,40:110)));
            %         Corrected_p_values(coherence,decoding_type,40:110)=mafdr(squeeze(p(coherence,decoding_type,40:110)));
        end
    end
    if region==1
        save('occipito_temporal_decoding_stim_resp_aligned_coh_levels_st.mat','accuracy_coherences_st','Corrected_p_values_st')     
    elseif region==2
        save('fronto_parietal_decoding_stim_resp_aligned_coh_levels_st.mat','accuracy_coherences_st','Corrected_p_values_st')        
    elseif region==3
        save('whole_brain_decoding_stim_resp_aligned_coh_levels_st.mat','accuracy_coherences_st','Corrected_p_values_st')
    end   
    plot(squeeze(nanmean(accuracy_coherences_st(:,:,1,40:110),2))')
    figure;
    plot(squeeze(nanmean(Corrected_p_values_st(:,1,40:110),2))')
else
    Corrected_p_values_rp=ones(4,6,200);
    for coherence=1:4
        for decoding_type=1:6
            Corrected_p_values_rp(coherence,decoding_type,91:161)=(squeeze(p(coherence,decoding_type,91:161)));
            %         Corrected_p_values(coherence,decoding_type,40:110)=mafdr(squeeze(p(coherence,decoding_type,40:110)));
        end
    end
    if region==1
        save('occipito_temporal_decoding_stim_resp_aligned_coh_levels_rp.mat','accuracy_coherences_rp','Corrected_p_values_rp')     
    elseif region==2
        save('fronto_parietal_decoding_stim_resp_aligned_coh_levels_rp.mat','accuracy_coherences_rp','Corrected_p_values_rp')        
    elseif region==3
        save('whole_brain_decoding_stim_resp_aligned_coh_levels_rp.mat','accuracy_coherences_rp','Corrected_p_values_rp')
    end
    plot(squeeze(nanmean(accuracy_coherences_rp(:,:,1,91:161),2))')
    figure;
    plot(squeeze(nanmean(Corrected_p_values_rp(:,1,91:161),2))')
end



