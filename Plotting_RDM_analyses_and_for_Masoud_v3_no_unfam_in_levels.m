
clc;
clear all;
% close all;
region=1;
stim_resp=2;
% for region=1:3
%     for stim_resp=1:2
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
                load(['st_aligned_partialRDM_IMG_SP_New_randoming_pls_NoPart_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'_Chris.mat']);
            else
                load(['rp_aligned_partialRDM_IMG_SP_New_randoming_pls_NoPart_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'_Chris.mat']);
            end
            
            True_correlations_Fam_Unfam=ParCorrelations_Fam_Unfam;
            Random_correlations_Fam_Unfam=ParCorrelations_Fam_Unfam_random;
            True_correlations_Fam_Levels=ParCorrelations_Fam_Levels;
            Random_correlations_Fam_Levels=ParCorrelations_Fam_Levels_random;
            
            %     Coh_Fam_Unfam(coherence)=plot(x,smooth(nanmean(True_correlations_Fam_Unfam),5),'Color',colors{coherence});
            %     hold on;
            Coh_Fam_Levels(coherence)=plot(x,smooth(nanmean(True_correlations_Fam_Levels),5),'LineStyle','--','Color',colors{coherence});
            hold on;
            for time=1:71
                if sum(time==inclusion_time)>0 && nanmean(True_correlations_Fam_Unfam(:,time),1)>0 && sum(nanmean(True_correlations_Fam_Unfam(:,time),1)>nanmean(Random_correlations_Fam_Unfam(:,time,:),1))>((1-p_value).*2000)
                    Significances_Fam_Unfam(coherence,time)=1;
                else
                    Significances_Fam_Unfam(coherence,time)=nan;
                end
                if sum(time==inclusion_time)>0 && nanmean(True_correlations_Fam_Levels(:,time),1)>0 && sum(nanmean(True_correlations_Fam_Levels(:,time),1)>nanmean(Random_correlations_Fam_Levels(:,time,:),1))>((1-p_value).*2000)
                    Significances_Fam_Levels(coherence,time)=1;
                else
                    Significances_Fam_Levels(coherence,time)=nan;
                end
            end
            %     plot(x,-0.002-Significances_Fam_Unfam(coherence,:).*coherence.*0.001,'*','Color',colors{coherence})
            %     plot(x,-0.006-Significances_Fam_Levels(coherence,:).*coherence.*0.001,'o','Color',colors{coherence})
            
            
            True_correlations_Fam_Unfam_NP=ParCorrelations_Fam_Unfam_NP;
            Random_correlations_Fam_Unfam_NP=ParCorrelations_Fam_Unfam_random_NP;
            True_correlations_Fam_Levels_NP=ParCorrelations_Fam_Levels_NP;
            Random_correlations_Fam_Levels_NP=ParCorrelations_Fam_Levels_random_NP;
            
            %     Coh_Fam_Unfam_NP(coherence)=plot(x,smooth(nanmean(True_correlations_Fam_Unfam_NP),5),'Color',colors{coherence});
            %     hold on;
            Coh_Fam_Levels_NP(coherence)=plot(x,smooth(nanmean(True_correlations_Fam_Levels_NP),5),'LineStyle','--','Color',colors{coherence});
            hold on;
            for time=1:71
                if sum(time==inclusion_time)>0 && nanmean(True_correlations_Fam_Unfam_NP(:,time),1)>0 && sum(nanmean(True_correlations_Fam_Unfam_NP(:,time),1)>nanmean(Random_correlations_Fam_Unfam_NP(:,time,:),1))>((1-p_value).*2000)
                    Significances_Fam_Unfam_NP(coherence,time)=1;
                else
                    Significances_Fam_Unfam_NP(coherence,time)=nan;
                end
                if sum(time==inclusion_time)>0 && nanmean(True_correlations_Fam_Levels_NP(:,time),1)>0 && sum(nanmean(True_correlations_Fam_Levels_NP(:,time),1)>nanmean(Random_correlations_Fam_Levels_NP(:,time,:),1))>((1-p_value).*2000)
                    Significances_Fam_Levels_NP(coherence,time)=1;
                else
                    Significances_Fam_Levels_NP(coherence,time)=nan;
                end
            end
            %     plot(x,-0.002-Significances_Fam_Unfam_NP(coherence,:).*coherence.*0.001,'*','Color',colors{coherence})
            %     plot(x,-0.006-Significances_Fam_Levels_NP(coherence,:).*coherence.*0.001,'o','Color',colors{coherence})
            
            
            %% saving for Masoud
            significance_Fam_Unfam=Significances_Fam_Unfam(coherence,:);
            significance_Fam_Levels=Significances_Fam_Levels(coherence,:);
            significance_Fam_Unfam_NP=Significances_Fam_Unfam_NP(coherence,:);
            significance_Fam_Levels_NP=Significances_Fam_Levels_NP(coherence,:);
            
            if stim_resp==1
                save(['st_aligned_partialRDM_New_summarized_for_Masoud_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'_Chris.mat'],'ParCorrelations_Fam_Unfam','ParCorrelations_Fam_Levels','significance_Fam_Unfam','significance_Fam_Levels','ParCorrelations_Fam_Unfam_NP','ParCorrelations_Fam_Levels_NP','significance_Fam_Unfam_NP','significance_Fam_Levels_NP');
            else
                save(['rp_aligned_partialRDM_New_summarized_for_Masoud_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'_Chris.mat'],'ParCorrelations_Fam_Unfam','ParCorrelations_Fam_Levels','significance_Fam_Unfam','significance_Fam_Levels','ParCorrelations_Fam_Unfam_NP','ParCorrelations_Fam_Levels_NP','significance_Fam_Unfam_NP','significance_Fam_Levels_NP');
            end
        end
%     end
% end