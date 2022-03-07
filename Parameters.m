clc;
clear all;
% close all;

%% Settings
% cats=0; % categories=1 and coh=0
only_answered_trials=0; % all trials=0
only_correct_trials=0; % plus incorrect=0 you should also change the saving file names
only_incorrect_trials=0; % plus correct=0

% stim_resp=1; % for stim=1, for resp=2
baseline_correction=0;
%% ERPs
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
sessions=[1 2];


for region=[3 2 1]
    if region==3
        Means=nan*ones(62,4,195,21);
        Vars=nan*ones(62,4,195,21);
        Cors=nan*ones(62,62,4,195,21);
    else
        Means=nan*ones(14,4,195,21);
        Vars=nan*ones(14,4,195,21);
        Cors=nan*ones(14,14,4,195,21);
    end
    for stim_resp=[1 2]
        for subj=[1:16 20 21]
            for coherence=[1:4]
                coherences=[0.22 0.3 0.45 0.55];
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
                        indx{cc,session}=([Exp_data{1,session}.stim.stimTrain.imageCategory]==counte & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
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

                %% Params
                time=0;
                steps=10;
                windows=50;
                for tws=1:steps:2000-windows
                    tw=tws:tws+windows;
                    time=time+1;
                    for cond=1:4
                        Means(:,cond,time,subj)=nanmean(squeeze(nanmean(signals{1,cond}(:,tw,:),2))');
                        Vars(:,cond,time,subj)=nanvar(squeeze(nanmean(signals{1,cond}(:,tw,:),2))');
                        for ch1=1:size(signals{1,cond},1)
                            for ch2=ch1+1:size(signals{1,cond},1)
                                Cors(ch1,ch2,cond,time,subj)=corr(squeeze(nanmean(signals{1,cond}(ch1,tw,:),2)),squeeze(nanmean(signals{1,cond}(ch2,tw,:),2)));
                            end
                        end
                    end
                    [region stim_resp subj coherence time]
                end

                if stim_resp==1
                    save(['st_ParametersV2_All_trials_region_',num2str(region),'_Coh',num2str(coherences(coherence)),'.mat'],'Means','Vars','Cors');
                else
                    save(['rp_ParametersV2_All_trials_region_',num2str(region),'_Coh',num2str(coherences(coherence)),'.mat'],'Means','Vars','Cors');
                end
                clearvars -except Means Vars Cors region test_on_errors sizes_signals accuracy Exp_data baseline_correction EEG_signals subject sessions steps baseline_span signals2 signals subj cats stim_resp only_answered_trials only_correct_trials only_incorrect_trials only_answered_trials2 only_correct_trials2 only_incorrect_trials2
            end           
        end
    end
end

