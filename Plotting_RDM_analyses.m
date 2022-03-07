
clc;
clear all;
% close all;
region=3;
stim_resp=2;

Coherences=[0.22 0.30 0.45 0.55];

if stim_resp==1
    x=[-100:10:600];
    inclusion_time=11:71;
else
    x=[-600:10:100];
    inclusion_time=1:71;
end
Significances=nan.*ones(4,71);
colors={'r','g','b','k'};
p_value=0.05;
load('RDM_matrices.mat');
figure;
for coherence=1:4
    if stim_resp==1
        load(['st_aligned_RDM_New_randoming_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'Correlations_Fam_Unfam','Correlations_Fam_Levels','Correlations_Fam_Unfam_random','Correlations_Fam_Levels_random');
    else
        load(['rp_aligned_RDM_New_randoming_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'Correlations_Fam_Unfam','Correlations_Fam_Levels','Correlations_Fam_Unfam_random','Correlations_Fam_Levels_random');
    end
    
    True_correlations_Fam_Unfam=Correlations_Fam_Unfam;
    Random_correlations_Fam_Unfam=Correlations_Fam_Unfam_random;
    True_correlations_Fam_Levels=Correlations_Fam_Levels;
    Random_correlations_Fam_Levels=Correlations_Fam_Levels_random;
    
    Coh_Fam_Unfam(coherence)=plot(x,nanmean(True_correlations_Fam_Unfam),'Color',colors{coherence});
    hold on;
    Coh_Fam_Levels(coherence)=plot(x,nanmean(True_correlations_Fam_Levels),'LineStyle','--','Color',colors{coherence});
    
    for time=1:71
        if sum(time==inclusion_time)>0 && nanmean(True_correlations_Fam_Unfam(:,time),1)>0 && sum(nanmean(True_correlations_Fam_Unfam(:,time),1)>nanmean(Random_correlations_Fam_Unfam(:,time,:),1))>((1-p_value).*10000)
            Significances_Fam_Unfam(coherence,time)=1;
        else
            Significances_Fam_Unfam(coherence,time)=nan;            
        end
        if sum(time==inclusion_time)>0 && nanmean(True_correlations_Fam_Levels(:,time),1)>0 && sum(nanmean(True_correlations_Fam_Levels(:,time),1)>nanmean(Random_correlations_Fam_Levels(:,time,:),1))>((1-p_value).*10000)
            Significances_Fam_Levels(coherence,time)=1;
        else
            Significances_Fam_Levels(coherence,time)=nan;
        end
    end
    plot(x,-0.002-Significances_Fam_Unfam(coherence,:).*coherence.*0.001,'*','Color',colors{coherence})
    plot(x,-0.006-Significances_Fam_Levels(coherence,:).*coherence.*0.001,'o','Color',colors{coherence})

    %% saving for Masoud
    significance_Fam_Unfam=Significances_Fam_Unfam(coherence,:);
    significance_Fam_Levels=Significances_Fam_Levels(coherence,:);
    if stim_resp==1
        save(['st_aligned_RDM_New_summarized_for_Masoud_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'Correlations_Fam_Unfam','Correlations_Fam_Levels','significance_Fam_Unfam','significance_Fam_Levels','Familiar_Unfamiliar_Model_RDM','Familiarity_levels_Model_RDM');
    else
        save(['rp_aligned_RDM_New_summarized_for_Masoud_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'Correlations_Fam_Unfam','Correlations_Fam_Levels','significance_Fam_Unfam','significance_Fam_Levels','Familiar_Unfamiliar_Model_RDM','Familiarity_levels_Model_RDM');
    end
end

