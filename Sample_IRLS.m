clc;
clear all;
close all;

%% Settings
% cats=0; % categories=1 and coh=0
only_answered_trials=1; % all trials=0
only_correct_trials=1; % plus incorrect=0
% stim_resp=1; % for stim=1, for resp=2
% subjs=[1 3 5];      % subj #
%% ERPs
for subj=[1:11]
    for cats=[0 1]
        for stim_resp=[1 2]
            load('All_subj_preprocessed_datasets_baseline_intact.mat');
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
            % if cats==1
            %     legend ('Control','Famous','Familiar','Self','Location','southeast');
            % else
            %     legend ('Noise = 0.22','Noise = 0.30','Noise = 0.45','Noise = 0.55','Location','southeast');
            % end
            %%
            % EOG AND M1 are channels #32 and #13
            for i=1:4
                signals{1,i}(13,:,:)=[];
                signals{1,i}(32,:,:)=[];
            end
            clearvars -except signals subj cats stim_resp only_answered_trials only_correct_trials
            %% IRLS
            landa=0.1;
            b=0.1;
            span=50;
            cond1=1;
            cond2=4;
            c=0;
            for tw=1:2000-span
                c=c+1;
                w=[zeros(size(signals{1,1},1),1);b];
                wind=tw:tw+span;
                wtb=nan.*ones(length(w),size(signals{1,cond1},3),size(signals{1,cond2},3));
                for tr1=1:size(signals{1,cond1},3)
                    for tr2=1:size(signals{1,cond2},3)
                        w=[zeros(size(signals{1,1},1),1);b];
                        X=[signals{1,cond1}(:,wind,tr1) signals{1,cond2}(:,wind,tr2)];
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
                %     wp(:,c)=w;
                wp(:,c)=squeeze(nanmean(nanmean(wtb,2),3));
                [subj cats stim_resp]
            end
            % for i=1:size(wp,2)-1
            %     diff(i)=abs(mean(wp(:,i+1)-wp(:,i)));
            % end
            % plot(diff)
            
            if stim_resp==1
                if cats==0
                    save(['st_al_dis_weights_coh_level_1_4_subj_',num2str(subj),'bs_in.mat'],'wp','span','landa','b');
                else
                    save(['st_al_dis_weights_cat_level_1_4_subj_',num2str(subj),'bs_in.mat'],'wp','span','landa','b');
                end
            else
                if cats==0
                    save(['rp_al_dis_weights_coh_level_1_4_subj_',num2str(subj),'bs_in.mat'],'wp','span','landa','b');
                else
                    save(['rp_al_dis_weights_cat_level_1_4_subj_',num2str(subj),'bs_in.mat'],'wp','span','landa','b');
                end
            end
            clearvars wp wt wtb
        end
    end
end









