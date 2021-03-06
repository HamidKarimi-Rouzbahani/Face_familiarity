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
iterations=1;
%% ERPs
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
sessions=[1 2];
ParCorrelations_Fam_Unfam_frnt=nan*ones(21,71);
ParCorrelations_Fam_Unfam_ocpt=nan*ones(21,71);
ParCorrelations_Fam_Levels_frnt=nan*ones(21,71);
ParCorrelations_Fam_Levels_ocpt=nan*ones(21,71);
ParCorrelations_FF=nan*ones(21,71);
ParCorrelations_FB=nan*ones(21,71);

ParCorrelations_FF_random=nan*ones(21,71,iterations);
ParCorrelations_FB_random=nan*ones(21,71,iterations);
ParCorrelations_Fam_Unfam_random_frnt=nan*ones(21,71,iterations);
ParCorrelations_Fam_Levels_random_frnt=nan*ones(21,71,iterations);
ParCorrelations_Fam_Unfam_random_ocpt=nan*ones(21,71,iterations);
ParCorrelations_Fam_Levels_random_ocpt=nan*ones(21,71,iterations);

Coherences=[0.22 0.30 0.45 0.55];
image_directory='D:\Hamid\Farzad\subject_images\Subject_';
for past_time=[170 200 400]
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
                signals_frnt=signals;
                signals_ocpt=signals;
                for i=1:4
                    signals_frnt{1,i}([13:32 44:57 60:64],:,:)=[];
                end
                for i=1:4
                    signals_ocpt{1,i}([1:19 32:47 58:59],:,:)=[];
                end
                
                % Adding all familiar categories to cond #2
                sizes_signals_frnt=[size(signals_frnt{1,1},3) size(signals_frnt{1,2},3) size(signals_frnt{1,3},3) size(signals_frnt{1,4},3)];
                signals_frnt{1,1}(:,:,end+1:end+size(signals_frnt{1,2},3))=signals_frnt{1,2};
                signals_frnt{1,1}(:,:,end+1:end+size(signals_frnt{1,3},3))=signals_frnt{1,3};
                signals_frnt{1,1}(:,:,end+1:end+size(signals_frnt{1,4},3))=signals_frnt{1,4};
                signals_frnt{1,2}=[];
                signals_frnt{1,3}=[];
                signals_frnt{1,4}=[];
                
                sizes_signals_ocpt=[size(signals_ocpt{1,1},3) size(signals_ocpt{1,2},3) size(signals_ocpt{1,3},3) size(signals_ocpt{1,4},3)];
                signals_ocpt{1,1}(:,:,end+1:end+size(signals_ocpt{1,2},3))=signals_ocpt{1,2};
                signals_ocpt{1,1}(:,:,end+1:end+size(signals_ocpt{1,3},3))=signals_ocpt{1,3};
                signals_ocpt{1,1}(:,:,end+1:end+size(signals_ocpt{1,4},3))=signals_ocpt{1,4};
                signals_ocpt{1,2}=[];
                signals_ocpt{1,3}=[];
                signals_ocpt{1,4}=[];
                clearvars -except iterations ParCorrelations_Fam_Unfam_frnt ParCorrelations_Fam_Unfam_ocpt ParCorrelations_Fam_Levels_frnt ParCorrelations_Fam_Levels_ocpt ParCorrelations_Fam_Unfam_random_frnt ParCorrelations_Fam_Levels_random_frnt ParCorrelations_Fam_Unfam_random_ocpt ParCorrelations_Fam_Levels_random_ocpt ParCorrelations_FF ParCorrelations_FB ParCorrelations_FF_random ParCorrelations_FB_random past_time past_window image_directory Coherences coherence region test_on_errors sizes_signals_frnt sizes_signals_ocpt accuracy Exp_data baseline_correction EEG_signals subject sessions steps baseline_span signals2 signals_frnt signals_ocpt subj cats stim_resp only_answered_trials only_correct_trials only_incorrect_trials only_answered_trials2 only_correct_trials2 only_incorrect_trials2
                %% Model RDMs
                Familiar_Unfamiliar_Model_RDM=nan*ones(size(signals_frnt{1,1},3),size(signals_frnt{1,1},3));
                Familiarity_levels_Model_RDM=nan*ones(size(signals_frnt{1,1},3),size(signals_frnt{1,1},3));
                for i=1:size(signals_frnt{1,1},3)
                    for j=i+1:size(signals_frnt{1,1},3)
                        if j<sizes_signals_frnt(1)+1 | i>sizes_signals_frnt(1)
                            Familiar_Unfamiliar_Model_RDM(i,j)=1;
                        else
                            Familiar_Unfamiliar_Model_RDM(i,j)=0;
                        end
                        
                        if j<sizes_signals_frnt(1)+1 | (i>sizes_signals_frnt(1) & j>sizes_signals_frnt(1) & j<(sizes_signals_frnt(1)+sizes_signals_frnt(2))+1) | (i>(sizes_signals_frnt(1)+sizes_signals_frnt(2)) & j>(sizes_signals_frnt(1)+sizes_signals_frnt(2)) & j<(sizes_signals_frnt(1)+sizes_signals_frnt(2)+sizes_signals_frnt(3)+1)) || (i>(sizes_signals_frnt(1)+sizes_signals_frnt(2)+sizes_signals_frnt(3)) & j>(sizes_signals_frnt(1)+sizes_signals_frnt(2)+sizes_signals_frnt(3)) & j<(sizes_signals_frnt(1)+sizes_signals_frnt(2)+sizes_signals_frnt(3)+sizes_signals_frnt(4)+1))
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
                %             past_time=125;
                past_window=50;
                
                if stim_resp==1
                    spans=400:steps:1100;
                else
                    spans=900:steps:1600;
                end
                for tws=spans
                    tw=tws:tws+steps;
                    time=time+1;
                    
                    %% Regular partial correlation for decoding
                    Neural_RDM_frnt_present=nan*ones(size(signals_frnt{1,1},3),size(signals_frnt{1,1},3));
                    for i=1:size(signals_frnt{1,1},3)
                        for j=i+1:size(signals_frnt{1,1},3)
                            Neural_RDM_frnt_present(i,j)=corr(squeeze(mean(signals_frnt{1,1}(:,tw,i),2)),squeeze(mean(signals_frnt{1,1}(:,tw,j),2)));
                        end
                    end
                    
                    Neural_RDM_ocpt_present=nan*ones(size(signals_ocpt{1,1},3),size(signals_ocpt{1,1},3));
                    for i=1:size(signals_ocpt{1,1},3)
                        for j=i+1:size(signals_ocpt{1,1},3)
                            Neural_RDM_ocpt_present(i,j)=corr(squeeze(mean(signals_ocpt{1,1}(:,tw,i),2)),squeeze(mean(signals_ocpt{1,1}(:,tw,j),2)));
                        end
                    end
                    
                    ParCorrelations_Fam_Unfam_frnt(subj,time)=partialcorr(reshape(Familiar_Unfamiliar_Model_RDM,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Neural_RDM_frnt_present,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Image_RDM,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),'rows','complete');
                    ParCorrelations_Fam_Levels_frnt(subj,time)=partialcorr(reshape(Familiarity_levels_Model_RDM,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Neural_RDM_frnt_present,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Image_RDM,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),'rows','complete');
                    
                    ParCorrelations_Fam_Unfam_ocpt(subj,time)=partialcorr(reshape(Familiar_Unfamiliar_Model_RDM,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Neural_RDM_ocpt_present,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Image_RDM,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),'rows','complete');
                    ParCorrelations_Fam_Levels_ocpt(subj,time)=partialcorr(reshape(Familiarity_levels_Model_RDM,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Neural_RDM_ocpt_present,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Image_RDM,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),'rows','complete');
                    
                    
                    inds_non_nans=find(~isnan(Familiar_Unfamiliar_Model_RDM));
                    
                    for it=1:iterations
                        Familiar_Unfamiliar_Model_RDM_tmp=nan*ones(size(signals_frnt{1,1},3));
                        Familiar_Unfamiliar_Model_RDM_tmp(randsample(inds_non_nans,sum(sum(Familiar_Unfamiliar_Model_RDM==1))))=1;
                        Familiar_Unfamiliar_Model_RDM_tmp([isnan(Familiar_Unfamiliar_Model_RDM_tmp) & ~isnan(Familiar_Unfamiliar_Model_RDM)])=0;
                        
                        Familiarity_levels_Model_RDM_tmp=nan*ones(size(signals_frnt{1,1},3));
                        Familiarity_levels_Model_RDM_tmp(randsample(inds_non_nans,sum(sum(Familiarity_levels_Model_RDM==1))))=1;
                        Familiarity_levels_Model_RDM_tmp([isnan(Familiarity_levels_Model_RDM_tmp) & ~isnan(Familiarity_levels_Model_RDM)])=0;
                        
                        ParCorrelations_Fam_Unfam_random_frnt(subj,time,it)=partialcorr(reshape(Familiar_Unfamiliar_Model_RDM_tmp,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Neural_RDM_frnt_present,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Image_RDM,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),'rows','complete');
                        ParCorrelations_Fam_Levels_random_frnt(subj,time,it)=partialcorr(reshape(Familiarity_levels_Model_RDM_tmp,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Neural_RDM_frnt_present,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Image_RDM,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),'rows','complete');
                        
                        ParCorrelations_Fam_Unfam_random_ocpt(subj,time,it)=partialcorr(reshape(Familiar_Unfamiliar_Model_RDM_tmp,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Neural_RDM_ocpt_present,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Image_RDM,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),'rows','complete');
                        ParCorrelations_Fam_Levels_random_ocpt(subj,time,it)=partialcorr(reshape(Familiarity_levels_Model_RDM_tmp,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Neural_RDM_ocpt_present,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Image_RDM,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),'rows','complete');
                    end
                    
                    %% Information flow analysis
                    
                    Neural_RDM_frnt_past=nan*ones(size(signals_frnt{1,1},3),size(signals_frnt{1,1},3));
                    for i=1:size(signals_frnt{1,1},3)
                        for j=i+1:size(signals_frnt{1,1},3)
                            Neural_RDM_frnt_past(i,j)=corr(squeeze(mean(signals_frnt{1,1}(:,[tw-past_time-0.5*past_window:tw-past_time+0.5*past_window],i),2)),squeeze(mean(signals_frnt{1,1}(:,[tw-past_time-0.5*past_window:tw-past_time+0.5*past_window],j),2)));
                        end
                    end
                    
                    Neural_RDM_ocpt_past=nan*ones(size(signals_ocpt{1,1},3),size(signals_ocpt{1,1},3));
                    for i=1:size(signals_ocpt{1,1},3)
                        for j=i+1:size(signals_ocpt{1,1},3)
                            Neural_RDM_ocpt_past(i,j)=corr(squeeze(mean(signals_ocpt{1,1}(:,[tw-past_time-0.5*past_window:tw-past_time+0.5*past_window],i),2)),squeeze(mean(signals_ocpt{1,1}(:,[tw-past_time-0.5*past_window:tw-past_time+0.5*past_window],j),2)));
                        end
                    end
                    
                    ParCorrelations_FF(subj,time)=partialcorr(reshape(Neural_RDM_frnt_present,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Neural_RDM_ocpt_past,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),[reshape(Neural_RDM_frnt_past,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]) reshape(Image_RDM,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1])],'rows','complete');
                    ParCorrelations_FB(subj,time)=partialcorr(reshape(Neural_RDM_ocpt_present,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Neural_RDM_frnt_past,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),[reshape(Neural_RDM_ocpt_past,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]) reshape(Image_RDM,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1])],'rows','complete');
                    
                    
                    inds_non_nans=find(~isnan(Neural_RDM_frnt_present));
                    
                    for it=1:iterations
                        Randomized_Neural_RDM_frnt_present=nan*ones(size(signals_frnt{1,1},3));
                        Randomized_Neural_RDM_frnt_present(randsample(inds_non_nans,length(inds_non_nans)))=Neural_RDM_frnt_present(inds_non_nans);
                        
                        Randomized_Neural_RDM_ocpt_present=nan*ones(size(signals_ocpt{1,1},3));
                        Randomized_Neural_RDM_ocpt_present(randsample(inds_non_nans,length(inds_non_nans)))=Neural_RDM_ocpt_present(inds_non_nans);
                        
                        ParCorrelations_FF_random(subj,time,it)=partialcorr(reshape(Randomized_Neural_RDM_frnt_present,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),reshape(Neural_RDM_ocpt_past,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),[reshape(Neural_RDM_frnt_past,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]) reshape(Image_RDM,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1])],'rows','complete');
                        ParCorrelations_FB_random(subj,time,it)=partialcorr(reshape(Randomized_Neural_RDM_ocpt_present,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]),reshape(Neural_RDM_frnt_past,[size(signals_frnt{1,1},3)*size(signals_frnt{1,1},3) 1]),[reshape(Neural_RDM_ocpt_past,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1]) reshape(Image_RDM,[size(signals_ocpt{1,1},3)*size(signals_ocpt{1,1},3) 1])],'rows','complete');
                    end
                    [past_time stim_resp coherence subj time]
                end
                
                if stim_resp==1
                    save(['st_al_pCor_IMG_occip_front_and_Flow_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam_frnt','ParCorrelations_Fam_Unfam_ocpt','ParCorrelations_Fam_Levels_frnt','ParCorrelations_Fam_Levels_ocpt','ParCorrelations_Fam_Unfam_random_frnt','ParCorrelations_Fam_Levels_random_frnt','ParCorrelations_Fam_Unfam_random_ocpt','ParCorrelations_Fam_Levels_random_ocpt','ParCorrelations_FF','ParCorrelations_FB','ParCorrelations_FF_random','ParCorrelations_FB_random');
                else
                    save(['rp_al_pCor_IMG_occip_front_and_Flow_pasttime_',num2str(past_time),'_pastwindow_',num2str(past_window),'_coherence_',num2str(Coherences(coherence)),'.mat'],'ParCorrelations_Fam_Unfam_frnt','ParCorrelations_Fam_Unfam_ocpt','ParCorrelations_Fam_Levels_frnt','ParCorrelations_Fam_Levels_ocpt','ParCorrelations_Fam_Unfam_random_frnt','ParCorrelations_Fam_Levels_random_frnt','ParCorrelations_Fam_Unfam_random_ocpt','ParCorrelations_Fam_Levels_random_ocpt','ParCorrelations_FF','ParCorrelations_FB','ParCorrelations_FF_random','ParCorrelations_FB_random');
                end
            end
        end
    end
end