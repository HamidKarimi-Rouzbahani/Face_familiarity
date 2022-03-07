clc;
clear all;
close all;

%% Settings
% cats=0; % categories=1 and coh=0
only_answered_trials=0; % all trials=0
only_correct_trials=0; % plus incorrect=0 you should also change the saving file names
only_incorrect_trials=0; % plus correct=0

only_answered_trials2=0; % all trials=0
only_correct_trials2=0; % plus incorrect=0 you should also change the saving file names
only_incorrect_trials2=0; % plus correct=0
% stim_resp=1; % for stim=1, for resp=2
baseline_correction=1;
steps=20;

%% ERPs
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
sessions=[1 2];
Correlations_Fam_Unfam=nan*ones(21,200);
Correlations_Fam_Levels=nan*ones(21,200);
Correlations_Fam_Unfam_random=nan*ones(21,200,100);
Correlations_Fam_Levels_random=nan*ones(21,200,100);
test_on_errors=0;
Coherences=[0.22 0.30 0.45 0.55];
for region=[1]
    for stim_resp=[1 2]
        for coherence=1:length(Coherences)
            for subj=[1:16 20 21]
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
                        colors=[0 0.45 0.7 1];
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
                
                
                %         if baseline_correction
                %             for i=1:4
                %                 signals{1,i}=signals{1,i}-repmat(nanmean(nanmean(signals{1,i}(:,baseline_span,:),3),2),[1 2000 size(signals{1,i},3)]);
                %             end
                %         end
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
                clearvars -except Correlations_Fam_Unfam_random Correlations_Fam_Levels_random Correlations_Fam_Levels Correlations_Fam_Unfam Coherences coherence region test_on_errors sizes_signals accuracy Exp_data baseline_correction EEG_signals subject sessions steps baseline_span signals2 signals subj cats stim_resp only_answered_trials only_correct_trials only_incorrect_trials only_answered_trials2 only_correct_trials2 only_incorrect_trials2
                
                %% model RDMs
                Familiar_Unfamiliar_Model_RDM=nan*ones(size(signals{1,1},3),size(signals{1,1},3));
                Familiarity_levels_Model_RDM=nan*ones(size(signals{1,1},3),size(signals{1,1},3));
                for i=1:size(signals{1,1},3)
                    for j=i:size(signals{1,1},3)
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
                %% Neural RDMs
                iterations=100;
                time=0;
                if stim_resp==1
                    spans=400:steps:1100-steps;
                else
                    spans=900:steps:1600-steps;
                end
                for tws=spans
                    tw=tws:tws+steps;
                    time=time+1;
                    Neural_RDM=nan*ones(size(signals{1,1},3),size(signals{1,1},3));
                    for i=1:size(signals{1,1},3)
                        for j=i:size(signals{1,1},3)
                            Neural_RDM(i,j)=corr(squeeze(mean(signals{1,1}(:,tw,i),2)),squeeze(mean(signals{1,1}(:,tw,j),2)));
                        end
                    end
                    Correlations_Fam_Unfam(subj,time)=corr(reshape(Familiar_Unfamiliar_Model_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Neural_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),'rows','complete');
                    Correlations_Fam_Levels(subj,time)=corr(reshape(Familiarity_levels_Model_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Neural_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),'rows','complete');
                    
                    
                    for it=1:iterations
                        rand_indices=randsample([1:size(signals{1,1},3)],size(signals{1,1},3));
                        Shuffled_signal=squeeze(mean(signals{1,1}(:,tw,rand_indices),2));
                        Random_neural_RDM=nan*ones(size(signals{1,1},3),size(signals{1,1},3));
                        for i=1:size(signals{1,1},3)
                            for j=i:size(signals{1,1},3)
                                Random_neural_RDM(i,j)=corr(Shuffled_signal(:,i),Shuffled_signal(:,j));
                            end
                        end
                        Correlations_Fam_Unfam_random(subj,time,it)=corr(reshape(Familiar_Unfamiliar_Model_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Random_neural_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),'rows','complete');
                        Correlations_Fam_Levels_random(subj,time,it)=corr(reshape(Familiarity_levels_Model_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),reshape(Random_neural_RDM,[size(signals{1,1},3)*size(signals{1,1},3) 1]),'rows','complete');
                    end
                    [region stim_resp coherence subj time]
                end
                
                if stim_resp==1
                    save(['st_aligned_RDM_corrrelations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'Correlations_Fam_Unfam','Correlations_Fam_Levels','Correlations_Fam_Unfam_random','Correlations_Fam_Levels_random');
                else
                    save(['rp_aligned_RDM_corrrelations_baselined_windowed_region_',num2str(region),'_coherence_',num2str(Coherences(coherence)),'.mat'],'Correlations_Fam_Unfam','Correlations_Fam_Levels','Correlations_Fam_Unfam_random','Correlations_Fam_Levels_random');
                end
            end
        end
    end
end
