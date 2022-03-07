clc;
clear all;
close all;

%% Settings
% cats=0; % categories=1 and coh=0
only_answered_trials=1; % all trials=0
only_correct_trials=0; % plus incorrect=0 you should also change the saving file names
only_incorrect_trials=0; % plus correct=0

% stim_resp=1; % for stim=1, for resp=2
baseline_correction=0;
steps=10;
famous_effect=0;
%% ERPs
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
sessions=[1 2];
for subj=[1:16 20 21]
    for cats=[1]
        for stim_resp=[1 2]
            if stim_resp==1
                baseline_span=[1:500];
            else
                baseline_span=[300:800];
            end
            load(['Subject_',subject{subj},'_preprosessed.mat']);
            cc=0;
            if cats==1
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
                        indx{cc,session}=[Exp_data{1,session}.stim.stimTrain.imageCategory]==counte;
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
                            EEG_signals{stim_resp,session}(31,:,:)=nan;
                        elseif subj==10 && session==1
                            EEG_signals{stim_resp,session}(7,:,:)=nan;
                        elseif subj==11 && session==1
                            EEG_signals{stim_resp,session}(41,:,:)=nan;
                            EEG_signals{stim_resp,session}(63,:,:)=nan;
                        end
                        tmpt=cat(3,tmpt,EEG_signals{stim_resp,session}(:,:,indx{cc,session}));
                    end
                    signals{cc}=tmpt;
                end
            else
                for counte=[0.22 0.3 0.45 0.55]
                    tmpt=[];
                    if subj==13
                        sessions=1;
                    else
                        sessions=[1 2];
                    end
                    cc=cc+1;
                    for session=sessions
                        colors=[0 0.45 0.7 1];
                        indx{cc,session}=[Exp_data{1,session}.stim.stimTrain.imageNoise]==counte;
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
                            EEG_signals{stim_resp,session}(31,:,:)=nan;
                        elseif subj==10 && session==1
                            EEG_signals{stim_resp,session}(7,:,:)=nan;
                        elseif subj==11 && session==1
                            EEG_signals{stim_resp,session}(41,:,:)=nan;
                            EEG_signals{stim_resp,session}(63,:,:)=nan;
                        end
                        tmpt=cat(3,tmpt,EEG_signals{stim_resp,session}(:,:,indx{cc,session}));
                    end
                    signals{cc}=tmpt;
                end
            end
            
            if baseline_correction
                for i=1:4
                    signals{1,i}=signals{1,i}-repmat(nanmean(nanmean(signals{1,i}(:,baseline_span,:),3),2),[1 2000 size(signals{1,i},3)]);
                end
            end
            % if cats==1
            %     legend ('Control','Famous','Familiar','Self','Location','southeast');
            % else
            %     legend ('Coherence = 0.22','Coherence = 0.30','Coherence = 0.45','Coherence = 0.55','Location','southeast');
            % end
            %%
            % EOG AND M1 are channels #32 and #13
            for i=1:4
                signals{1,i}(13,:,:)=[];
                signals{1,i}(32,:,:)=[];
            end
            %% Adding all familiar categories to cond #2
            signals{1,2}(:,:,end+1:end+size(signals{1,3},3))=signals{1,3};
            signals{1,2}(:,:,end+1:end+size(signals{1,4},3))=signals{1,4};
            clearvars -except famous_effect subject sessions steps baseline_correction baseline_span signals subj cats stim_resp only_answered_trials only_correct_trials only_incorrect_trials
            %% IRLS
            landa=0.1;
            b=0.1;
            span=50;
            if famous_effect==0
                cond1=1;
                cond2=2;
            end
            c=0;
            for tw=1:steps:2000-span
                c=c+1;
                w=[zeros(size(signals{1,1},1),1);b];
                wind=tw:tw+span;
                wtb=nan.*ones(length(w),size(signals{1,cond1},3),size(signals{1,cond2},3));
                for tr1=1:size(signals{1,cond1},3)
                    for tr2=1:size(signals{1,cond2},3)
                        w=[zeros(size(signals{1,1},1),1);b];
                        X=[signals{1,cond1}(:,wind,tr1) signals{1,cond2}(:,wind,tr2)];
                        X(isnan(X))=0;
                        d=[ones([length(wind),1]);zeros([length(wind),1])];
                        Landa=landa.*[eye(size(X,1)) zeros(size(X,1),1);zeros(1,size(X,1)+1)];
                        X=[X;ones(1,size(X,2))];
                        p=1./(1+exp(-w'*X))';
                        g=X*(d-p)-Landa*w;
                        H=X*diag(p.*(1-p))*X'+Landa*eye(size(X,1));
                        w=w+pinv(H)*g;
                        wtb(:,tr1,tr2)=w;
                    end
                end
                wp(:,c)=squeeze(nanmean(nanmean(wtb,2),3));
                [subj cats stim_resp tw]
            end
            
            if stim_resp==1
                    save(['st_al_weights_fam_unf_subj_',num2str(subj),'conc_Alltrials.mat'],'wp','span','landa','b');
            else
                    save(['rp_al_weights_fam_unf_subj_',num2str(subj),'conc_Alltrials.mat'],'wp','span','landa','b');
            end
            clearvars wp wt wtb
        end
    end
end









