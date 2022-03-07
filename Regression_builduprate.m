clc;
clear all;
close all;
%% Loading and analysis of weights
only_answered_trials=1; % all trials=0
only_correct_trials=1; % plus incorrect=0
for subj=[1]
    for cats=[0]
        for stim_resp=[1]
            load('All_subj_preprocessed_datasets.mat');
            cc=0;
            if cats==1
                for counte=[1 2 4 3]
                    colors=[0 0.45 0.7 1];
                    cc=cc+1;
                    indx{cc,subj}=[Exp_data{1,subj}.stim.stimTrain.imageCategory]==counte;
                    if only_answered_trials==1
                        tmp=indx{cc,subj};
                        tmp(isnan(Exp_data{1,subj}.stim.ResponseData.Values(2,:)))=0;
                        indx{cc,subj}=tmp;
                        if only_correct_trials==1
                            tmp=indx{cc,subj};
                            tmp(Exp_data{1,subj}.stim.ResponseData.Values(2,:)==0)=0;
                            indx{cc,subj}=tmp;
                        end
                    end
                    signals{cc}=EEG_signals{stim_resp,subj}(:,:,indx{cc,subj});
                end
            else
                for counte=[0.22 0.3 0.45 0.55]
                    cc=cc+1;
                    indx{cc,subj}=[Exp_data{1,subj}.stim.stimTrain.imageNoise]==counte;
                    if only_answered_trials==1
                        tmp=indx{cc,subj};
                        tmp(isnan(Exp_data{1,subj}.stim.ResponseData.Values(2,:)))=0;
                        indx{cc,subj}=tmp;
                        if only_correct_trials==1
                            tmp=indx{cc,subj};
                            tmp(Exp_data{1,subj}.stim.ResponseData.Values(2,:)==0)=0;
                            indx{cc,subj}=tmp;
                        end
                    end
                    signals{cc}=EEG_signals{stim_resp,subj}(:,:,indx{cc,subj});
                end
            end
            %% EOG AND M1 are channels #32 and #13
            for i=1:4
                signals{1,i}(13,:,:)=[];
                signals{1,i}(32,:,:)=[];
            end
            clearvars -except signals subj cats stim_resp only_answered_trials only_correct_trials
            %%
            if stim_resp==1
                if cats==0
                    load(['st_al_dis_weights_coh_level_1_4_subj_',num2str(subj),'.mat'],'wp','span','landa','b');
                else
                    load(['st_al_dis_weights_cat_level_1_4_subj_',num2str(subj),'.mat'],'wp','span','landa','b');
                end
            else
                if cats==0
                    load(['rp_al_dis_weights_coh_level_1_4_subj_',num2str(subj),'.mat'],'wp','span','landa','b');
                else
                    load(['rp_al_dis_weights_cat_level_1_4_subj_',num2str(subj),'.mat'],'wp','span','landa','b');
                end
            end
            cond1=1;
            cond2=4;
            span=50;
            if stim_resp==1
                tws=1:1250-span;
            else
                tws=451:1700-span;
            end
            
            clearvars X
            c=0;
            d=[ones(1,size(signals{1,cond1}(:,1,:),3)) zeros(1,size(signals{1,cond2}(:,1,:),3))];
            for t=tws
                c=c+1;
                X=[mean(signals{1,cond1}(:,t,:),3) mean(signals{1,cond2}(:,t,:),3); 1 1];
                Xp=wp(:,c)'*squeeze(X);
                XX=[squeeze(signals{1,cond1}(:,t,:)) squeeze(signals{1,cond2}(:,t,:));ones(1,size(signals{1,cond1}(:,t,:),3)+size(signals{1,cond2}(:,t,:),3))];
                XXp=wp(:,c)'*XX;
                for trl=1:size(XXp,2)
                    for cls=1:2
                        dis(cls,trl)=abs(Xp(cls)-XXp(trl));
                    end
                    dis(:,trl)=dis(:,trl)./sum(dis(:,trl));
                end
                [~,~,~,AUC(c)] = perfcurve(d',1-dis(1,:)','1');
            end
            [~,indx_max_AUC]=max(AUC);
            
            ccc
            clearvars wp wt wtb
        end
    end
end