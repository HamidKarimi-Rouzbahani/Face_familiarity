%% Plotting
clc;
clear all;
close all;

stim_resp=1;
region=3;
% coherence=4;
Coherences=[0.22 0.30 0.45 0.55];



for coherence=[1:4]
    if stim_resp==1
        load(['st_aligned_partialRDM_IMG_New_randoming_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam','ParCorrelations_Fam_Levels','ParCorrelations_Fam_Unfam_random','ParCorrelations_Fam_Levels_random');
   else
        load(['rp_aligned_partialRDM_IMG_New_randoming_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam','ParCorrelations_Fam_Levels','ParCorrelations_Fam_Unfam_random','ParCorrelations_Fam_Levels_random');
   end
       plot(smooth(nanmean(ParCorrelations_Fam_Unfam),3));
       hold on;  
end
figure;
for coherence=[1:4]
    if stim_resp==1
        load(['st_aligned_partialRDM_IMG_New_randoming_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam','ParCorrelations_Fam_Levels','ParCorrelations_Fam_Unfam_random','ParCorrelations_Fam_Levels_random');
    else
        load(['rp_aligned_partialRDM_IMG_New_randoming_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam','ParCorrelations_Fam_Levels','ParCorrelations_Fam_Unfam_random','ParCorrelations_Fam_Levels_random');
    end
       plot(smooth(nanmean(ParCorrelations_Fam_Levels),3));
       hold on;    
end


%% Plotting
clc;
clear all;
close all;

stim_resp=2;
region=3;
Coherences=[0.22 0.30 0.45 0.55];
colors={[1 0 0],[0 1 0],[0 0 1],[0 0 0]};

for coherence=[1:4]
    if stim_resp==1
        load(['st_aligned_partialRDM_IMG_SP_New_randoming_pls_NoPart_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat']);
    else
        load(['rp_aligned_partialRDM_IMG_SP_New_randoming_pls_NoPart_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat']);
   end
       plot(smooth(nanmean(ParCorrelations_Fam_Unfam),10),'color',colors{coherence},'linestyle','--');
       hold on;
       plot(smooth(nanmean(ParCorrelations_Fam_Unfam_NP),10),'color',colors{coherence});       
end
figure;
for coherence=[1:4]
    if stim_resp==1
        load(['st_aligned_partialRDM_IMG_New_randoming_pls_NoPart_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat']);
    else
        load(['rp_aligned_partialRDM_IMG_New_randoming_pls_NoPart_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat']);
    end
       plot(smooth(nanmean(ParCorrelations_Fam_Levels),10),'color',colors{coherence});
       hold on;    
       plot(smooth(nanmean(ParCorrelations_Fam_Levels_NP),10),'color',colors{coherence},'linestyle','--');       
end


