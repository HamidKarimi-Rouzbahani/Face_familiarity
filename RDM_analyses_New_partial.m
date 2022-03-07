clc;
clear all;
close all;

%% Settings
% cats=0; % categories=1 and coh=0
only_answered_trials=1; % all trials=0
only_correct_trials=1; % plus incorrect=0 you should also change the saving file names
only_incorrect_trials=0; % plus correct=0

% stim_resp=1; % for stim=1, for resp=2
baseline_correction=1;
steps=10;

%% ERPs
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
sessions=[1 2];
iterations=10000;
ParCorrelations_Fam_Unfam=nan*ones(21,71);
ParCorrelations_Fam_Levels=nan*ones(21,71);
ParCorrelations_Fam_Unfam_random=nan*ones(21,71,iterations);
ParCorrelations_Fam_Levels_random=nan*ones(21,71,iterations);
test_on_errors=0;
Coherences=[0.22 0.30 0.45 0.55];
image_directory='C:\Users\mq20185770\Documents\MATLAB\Farzad\subject_images\Subject_';
for region=[3]
    for stim_resp=[1]
        for coherence=1:length(Coherences)
            for subj=[1:16]
                %% First dataset
                if stim_resp==1
                    baseline_span=[1:500];
                else
                    baseline_span=[300:800];
                end
                load(['Subject_',subject{subj},'_preprosessed.mat']);
                cc=0;
                for counte=[1 2 4 3]
                    tmpt=[];
                    cc=cc+1;
                    if subj==13
                        sessions=1;
                    else
                        sessions=[1 2];
                    end
                    for session=sessions
                        %% Filling in empty trial values
                        for i=1:size([Exp_data{1, session}.stim.stimTrain],2)
                            if isempty([Exp_data{1, session}.stim.stimTrain(i).imageCategory])==1
                                Exp_data{1, session}.stim.stimTrain(i).imageNoise=0;
                                Exp_data{1, session}.stim.stimTrain(i).imageCategory=0;
                            end
                        end
                        indx{cc,session}=([Exp_data{1,session}.stim.stimTrain.imageCategory]==counte & [Exp_data{1,session}.stim.stimTrain.imageNoise]==Coherences(coherence));
                        if only_answered_trials==1
                            tmp=indx{cc,session};
                            tmp(isnan(Exp_data{1,session}.stim.ResponseData.Values(2,:)))=0;
                            indx{cc,session}=tmp;
                            if only_correct_trials==1
                                tmp=indx{cc,session};
                                tmp(Exp_data{1,session}.stim.ResponseData.Values(2,:)==0)=0;
                                indx{cc,session}=tmp;
                            end
                            if only_incorrect_trials==1
                                tmp=indx{cc,session};
                                tmp(Exp_data{1,session}.stim.ResponseData.Values(2,:)==1)=0;
                                indx{cc,session}=tmp;
                            end
                        end
                        if subj==2 && session==2    % remove bad channels
                            EEG_signals{stim_resp,session}(31,:,:)=0;
                        elseif subj==10 && session==1
                            EEG_signals{stim_resp,session}(7,:,:)=0;
                        elseif subj==11 && session==1
                            EEG_signals{stim_resp,session}(41,:,:)=0;
                            EEG_signals{stim_resp,session}(63,:,:)=0;
                        end
                        tmpt=cat(3,tmpt,EEG_signals{stim_resp,session}(:,:,indx{cc,session}));
                    end
                    signals{cc}=tmpt;
                end
                
                
                if baseline_correction
                    tt=nanmean([nanmean(nanmean(signals{1,1}(:,baseline_span,:),3),2) nanmean(nanmean(signals{1,2}(:,baseline_span,:),3),2) nanmean(nanmean(signals{1,3}(:,baseline_span,:),3),2) nanmean(nanmean(signals{1,4}(:,baseline_span,:),3),2)],2);
                    for i=1:4
                        signals{1,i}=signals{1,i}-repmat(tt,[1 2000 size(signals{1,i},3)]);
                    end
                end
                
                % if cats==1
                %     legend ('Control','Famous','Familiar','Self','Location','southeast');
                % else
                %     legend ('Coherence = 0.22','Coherence = 0.30','Coherence = 0.45','Coherence = 0.55','Location','southeast');
                % end
                
                % EOG AND M1 are channels #32 and #13 plus non-region remove
                if region==1
                    for i=1:4
                        signals{1,i}([1:12 13 14:23 25:27 32 33:49 51:52 58:61],:,:)=[];
                    end
                elseif region==2
                    for i=1:4
                        signals{1,i}([1:5 7:9 12 13 14:15 17:20 23:25 27:31 32 33:37 40:41 43:44 47:50 53:64],:,:)=[];
                    end
                elseif region==3
                    for i=1:4
                        signals{1,i}([13 32],:,:)=[];
                    end
                end
                % Adding all familiar categories to cond #2
                sizes_signals=[size(signals{1,1},3) size(signals{1,2},3) size(signals{1,3},3) size(signals{1,4},3)];
                signals{1,1}(:,:,end+1:end+size(signals{1,2},3))=signals{1,2};
                signals{1,1}(:,:,end+1:end+size(signals{1,3},3))=signals{1,3};
                signals{1,1}(:,:,end+1:end+size(signals{1,4},3))=signals{1,4};
                signals{1,2}=[];
                signals{1,3}=[];
                signals{1,4}=[];
                clearvars -except iterations image_directory ParCorrelations_Fam_Unfam_random ParCorrelations_Fam_Levels_random ParCorrelations_Fam_Levels ParCorrelations_Fam_Unfam Coherences coherence region test_on_errors sizes_signals accuracy Exp_data baseline_correction EEG_signals subject sessions steps baseline_span signals2 signals subj cats stim_resp only_answered_trials only_correct_trials only_incorrect_trials only_answered_trials2 only_correct_trials2 only_incorrect_trials2
                %% Model RDMs
                Familiar_Unfamiliar_Model_RDM=nan*ones(size(signals{1,1},3),size(signals{1,1},3));
                Familiarity_levels_Model_RDM=nan*ones(size(signals{1,1},3),size(signals{1,1},3));
                for i=1:size(signals{1,1},3)
                    for j=i+1:size(signals{1,1},3)
                        if j<sizes_signals(1)+1 | i>sizes_signals(1)
                            Familiar_Unfamiliar_Model_RDM(i,j)=1;
                        else
                            Familiar_Unfamiliar_Model_RDM(i,j)=0;
                        end
                        
                        if j<sizes_signals(1)+1 | (i>sizes_signals(1) & j>sizes_signals(1) & j<(sizes_signals(1)+sizes_signals(2))+1) | (i>(sizes_signals(1)+sizes_signals(2)) & j>(sizes_signals(1)+sizes_signals(2)) & j<(sizes_signals(1)+sizes_signals(2)+sizes_signals(3)+1)) | (i>(sizes_signals(1)+sizes_signals(2)+sizes_signals(3)) & j>(sizes_signals(1)+sizes_signals(2)+sizes_signals(3)) & j<(sizes_signals(1)+sizes_signals(2)+sizes_signals(3)+sizes_signals(4)+1))
                            Familiarity_levels_Model_RDM(i,j)=1;
                        else
                            Familiarity_levels_Model_RDM(i,j)=0;
                        end
                    end
                end
                %% Image RDMs
                cc=0;
                images=[];
                for counte=[1 2 4 3]
                    tmpt=[];
                    cc=cc+1;
                    if subj==13
                        sessions=1;
                    else
                        sessions=[1 2];
                    end
                    for session=sessions
                        %% Filling in empty trial values
                        for i=1:size([Exp_data{1, session}.stim.stimTrain],2)
                            if isempty([Exp_data{1, session}.stim.stimTrain(i).imageCategory])==1
                                Exp_data{1, session}.stim.stimTrain(i).imageNoise=0;
                                Exp_data{1, session}.stim.stimTrain(i).imageCategory=0;
                            end
                        end
                        indx{cc,session}=([Exp_data{1,session}.stim.stimTrain.imageCategory]==counte & [Exp_data{1,session}.stim.stimTrain.imageNoise]==Coherences(coherence));
                        if only_answered_trials==1
                            tmp=indx{cc,session};
                            tmp(isnan(Exp_data{1,session}.stim.ResponseData.Values(2,:)))=0;
                            indx{cc,session}=tmp;
                            if only_correct_trials==1
                                tmp=indx{cc,session};
                                tmp(Exp_data{1,session}.stim.ResponseData.Values(2,:)==0)=0;
                                indx{cc,session}=tmp;
                            end
                            if only_incorrect_trials==1
                                tmp=indx{cc,session};
                                tmp(Exp_data{1,session}.stim.ResponseData.Values(2,:)==1)=0;
                                indx{cc,session}=tmp;
                            end
                        end
                        indices=find(indx{cc,session}==1);
                        for img=1:sum(indx{cc,session}==1)
                            loaded_img=reshape(imread([image_directory,subject{subj},'\coherence_',num2str(Coherences(coherence)*100),'\',Exp_data{1,session}.stim.ResponseData.Names{1,indices(img)}]),[400*400 1]);
                            tmpt=cat(2,tmpt,loaded_img);
                        end                        
                    end
                    images=cat(2,images,tmpt);
                end
                images=double(images);
                Image_RDM=nan*ones(size(images,2));               
                for i=1:size(images,2)
                    for j=i+1:size(images,2)
                        Image_RDM(i,j)=corr(images(:,i),images(:,j));
                    end
                end
                clearvars images tmpt
                %% Neural RDMs
                time=0;
                if stim_resp==1
                    spans=400:steps:1100;
                else
                    spans=900:steps:1600;
                end
                for tws=spans
                    tw=tws:tws+steps;
                    time=time+1;
                    Neural_RDM=nan*ones(size(signals{1,1},3),size(signals{1,1},3));
                    for i=1:size(signals{1,1},3)
                        for j=i+1:size(signals{1,1},3)
                            Neural_RDM(i,j)=corr(squeeze(mean(signals{1,1}(:,tw,i),2)),squeeze(mean(signals{1,1}(:,tw,j),2)));
                        end
                    end
                    ParCorrelations_Fam_Unfam(subj,time)=partialcorr(reshape(Familiar_Unfamiliar_Model_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Neural_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Image_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),'rows','complete');
                    ParCorrelations_Fam_Levels(subj,time)=partialcorr(reshape(Familiarity_levels_Model_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Neural_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Image_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),'rows','complete');
                    
                    inds_non_nans_FamUnfam=find(~isnan(Familiar_Unfamiliar_Model_RDM));
                    inds_non_nans_FamLev=find(~isnan(Familiarity_levels_Model_RDM));
                    
                    for it=1:iterations
                        Familiar_Unfamiliar_Model_RDM_tmp=nan*ones(size(signals{1,1},3),size(signals{1,1},3));
                        Familiar_Unfamiliar_Model_RDM_tmp(randsample(inds_non_nans_FamUnfam,sum(sum(Familiar_Unfamiliar_Model_RDM==1))))=1;
                        Familiar_Unfamiliar_Model_RDM_tmp([isnan(Familiar_Unfamiliar_Model_RDM_tmp) & ~isnan(Familiar_Unfamiliar_Model_RDM)])=0;
                        
                        Familiarity_levels_Model_RDM_tmp=nan*ones(size(signals{1,1},3),size(signals{1,1},3));
                        Familiarity_levels_Model_RDM_tmp(randsample(inds_non_nans_FamLev,sum(sum(Familiarity_levels_Model_RDM==1))))=1;
                        Familiarity_levels_Model_RDM_tmp([isnan(Familiarity_levels_Model_RDM_tmp) & ~isnan(Familiarity_levels_Model_RDM)])=0;
                        
                        ParCorrelations_Fam_Unfam_random(subj,time,it)=partialcorr(reshape(Familiar_Unfamiliar_Model_RDM_tmp,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Neural_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Image_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),'rows','complete');
                        ParCorrelations_Fam_Levels_random(subj,time,it)=partialcorr(reshape(Familiarity_levels_Model_RDM_tmp,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Neural_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Image_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),'rows','complete');
                    end
                    [region stim_resp coherence subj time]
                end
                ccc
                if stim_resp==1
                    save(['st_aligned_partialRDM_IMG_New_randoming_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam','ParCorrelations_Fam_Levels','ParCorrelations_Fam_Unfam_random','ParCorrelations_Fam_Levels_random');
                else
                    save(['rp_aligned_partialRDM_IMG_New_randoming_correlations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam','ParCorrelations_Fam_Levels','ParCorrelations_Fam_Unfam_random','ParCorrelations_Fam_Levels_random');
                end
                
            end
        end
    end
end
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



