clc;
clear all;
close all;
%% Loading and analysis of weights
only_answered_trials=1; % all trials=0
only_correct_trials=0; % plus incorrect=0
only_incorrect_trials=1; % plus correct=0 you should also change the loading file names

cats=1;
stim_resp=2;
famous_effect=0;

subjects=[1:16 20 21];
Xp_sorted_t=nan*ones(480,1950,length(subjects));
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
sessions=[1 2];
subs=0;
g=0;
for subj=subjects
    load(['Subject_',subject{subj},'_preprosessed.mat']);
    subs=subs+1;
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
    else
        for counte=[0.22 0.3 0.45 0.55]
            tmpt=[];
            cc=cc+1;
            if subj==13
                sessions=1;
            else
                sessions=[1 2];
            end
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
    end
    %% EOG AND M1 are channels #32 and #13
    for i=1:4
        signals{1,i}(13,:,:)=[];
        signals{1,i}(32,:,:)=[];
    end
    clearvars -except g subject subjects famous_effect sessions Exp_data indx_max_AUC Xp_sorted XXX XXXp AUC AUCn AUC_after AUCn_after a subs EEG_signals signals subj cats stim_resp only_answered_trials only_correct_trials only_incorrect_trials
    %% IRLS
    if famous_effect==0
        if only_correct_trials==1
            if stim_resp==1
                if cats==0
                    load(['st_al_weights_coh_level_1_4_subj_',num2str(subj),'conc.mat'],'wp','span','landa','b');
                else
                    load(['st_al_weights_cat_level_1_4_subj_',num2str(subj),'conc.mat'],'wp','span','landa','b');
                end
            else
                if cats==0
                    load(['rp_al_weights_coh_level_1_4_subj_',num2str(subj),'conc.mat'],'wp','span','landa','b');
                else
                    load(['rp_al_weights_cat_level_1_4_subj_',num2str(subj),'conc.mat'],'wp','span','landa','b');
                end
            end
        else
            if stim_resp==1
                if cats==0
                    load(['st_al_weights_coh_level_1_4_subj_',num2str(subj),'conc_ALLtrials.mat'],'wp','span','landa','b');
                else
                    load(['st_al_weights_cat_level_1_4_subj_',num2str(subj),'conc_ALLtrials.mat'],'wp','span','landa','b');
                end
            else
                if cats==0
                    load(['rp_al_weights_coh_level_1_4_subj_',num2str(subj),'conc_ALLtrials.mat'],'wp','span','landa','b');
                else
                    load(['rp_al_weights_cat_level_1_4_subj_',num2str(subj),'conc_ALLtrials.mat'],'wp','span','landa','b');
                end
            end
        end
    else
        if stim_resp==1
            if cats==0
                load(['st_al_weights_coh_level_1_4_subj_',num2str(subj),'conc_famous.mat'],'wp','span','landa','b');
            else
                load(['st_al_weights_cat_level_1_4_subj_',num2str(subj),'conc_famous.mat'],'wp','span','landa','b');
            end
        else
            if cats==0
                load(['rp_al_weights_coh_level_1_4_subj_',num2str(subj),'conc_famous.mat'],'wp','span','landa','b');
            else
                load(['rp_al_weights_cat_level_1_4_subj_',num2str(subj),'conc_famous.mat'],'wp','span','landa','b');
            end
        end
    end

    cond1=1;
    if famous_effect==1
        cond2=2;
    else
        cond2=4;
    end
    span=50;
    steps=10;
    tws=1:steps:2000-span;
    clearvars X
    if size(signals{1,cond2}(:,1,:),3)==0
        g=g+1;
        AUC(subs,:)=nan*ones(1,length(tws));
        AUCn(subs,:)=nan*ones(1,length(tws));
        AUC_after(subs,:)=nan*ones(1,length(tws));
        AUCn_after(subs,:)=nan*ones(1,length(tws));
        if stim_resp==1
            indx_max_AUC(subs)=95;
        else
            indx_max_AUC(subs)=140;
        end        
    else  
        c=0;
        d=[ones(1,size(signals{1,cond1}(:,1,:),3)) zeros(1,size(signals{1,cond2}(:,1,:),3))];
    for t=tws
        c=c+1;
        X=[mean(signals{1,cond1}(:,t,:),3) mean(signals{1,cond2}(:,t,:),3); 1 1];
        Xp=wp(:,c)'*squeeze(X);
        Xpn=mean(X);
        XX=[squeeze(signals{1,cond1}(:,t,:)) squeeze(signals{1,cond2}(:,t,:));ones(1,size(signals{1,cond1}(:,t,:),3)+size(signals{1,cond2}(:,t,:),3))];
        XXp=wp(:,c)'*XX;
        XXpn=mean(XX);
        for trl=1:size(XXp,2)
            for cls=1:2
                dis(cls,trl)=abs(Xp(cls)-XXp(trl));
                disn(cls,trl)=abs(Xpn(cls)-XXpn(trl));
            end
            dis(:,trl)=dis(:,trl)./sum(dis(:,trl));
            disn(:,trl)=disn(:,trl)./sum(disn(:,trl));
        end
        [~,~,~,AUC(subs,c)] = perfcurve(d',1-dis(1,:)','1');
        [~,~,~,AUCn(subs,c)] = perfcurve(d',1-disn(1,:)','1');
    end
    
    if famous_effect==0
        if stim_resp==1
            [~,indx_max_AUC(subs)]=max(AUC(subs,[850:steps:1050]./steps));
            indx_max_AUC(subs)=indx_max_AUC(subs)+850./steps;
%             [~,indx_max_AUC(subs)]=max(AUC(subs,[500:steps:1500]./steps));
%             indx_max_AUC(subs)=indx_max_AUC(subs)+500./steps;
         else
            [~,indx_max_AUC(subs)]=max(AUC(subs,[1300:steps:1500]./steps));
            indx_max_AUC(subs)=indx_max_AUC(subs)+1300./steps;
%             [~,indx_max_AUC(subs)]=max(AUC(subs,[500:steps:1700]./steps));
%             indx_max_AUC(subs)=indx_max_AUC(subs)+500./steps;
        end
    else
        if stim_resp==1
            [~,indx_max_AUC(subs)]=max(AUC(subs,[500:steps:1500]./steps));
            indx_max_AUC(subs)=indx_max_AUC(subs)+500./steps;
        else
            [~,indx_max_AUC(subs)]=max(AUC(subs,[500:steps:1700]./steps));
            indx_max_AUC(subs)=indx_max_AUC(subs)+500./steps;
        end
    end
    
    % AUC after applying best AUC time
    clearvars X
    c=0;
    for t=tws
        c=c+1;
        X=[mean(signals{1,cond1}(:,t,:),3) mean(signals{1,cond2}(:,t,:),3);1 1];
        Xp=wp(:,indx_max_AUC(subs))'*squeeze(X);
        Xpn=mean(X);
        XX=[squeeze(signals{1,cond1}(:,t,:)) squeeze(signals{1,cond2}(:,t,:));ones(1,size(signals{1,cond1}(:,t,:),3)+size(signals{1,cond2}(:,t,:),3))];
        XXp=wp(:,indx_max_AUC(subs))'*XX;
        XXpn=mean(XX);
        for trl=1:size(XXp,2)
            for cls=1:2
                dis(cls,trl)=abs(Xp(cls)-XXp(trl));
                disn(cls,trl)=abs(Xpn(cls)-XXpn(trl));
            end
            dis(:,trl)=dis(:,trl)./sum(dis(:,trl));
            disn(:,trl)=disn(:,trl)./sum(disn(:,trl));
        end
        [~,~,~,AUC_after(subs,c)] = perfcurve(d',1-dis(1,:)','1');
        [~,~,~,AUCn_after(subs,c)] = perfcurve(d',1-disn(1,:)','1');
    end
    end
    %% AUC significance
    %     iteration=10;
    %     for iter=1:iteration
    %         d=d(randperm(length(d)));
    %         c=0;
    %         for t=tws
    %             c=c+1;
    %             X=[mean(signals{1,cond1}(:,t,:),3) mean(signals{1,cond2}(:,t,:),3); 1 1];
    %             Xp=wp(:,c)'*squeeze(X);
    %             XX=[squeeze(signals{1,cond1}(:,t,:)) squeeze(signals{1,cond2}(:,t,:));ones(1,size(signals{1,cond1}(:,t,:),3)+size(signals{1,cond2}(:,t,:),3))];
    %             XXp=wp(:,c)'*XX;
    %             for trl=1:size(XXp,2)
    %                 for cls=1:2
    %                     dis(cls,trl)=abs(Xp(cls)-XXp(trl));
    %                 end
    %                 dis(:,trl)=dis(:,trl)./sum(dis(:,trl));
    %             end
    %             [~,~,~,AUC(c)] = perfcurve(d',1-dis(1,:)','1');
    %         end
    %         AUCit(iter,:)=AUC;
    %         iter
    %     end
    %
    %     sig_threshold=0.95;
    %     for i=1:size(AUCit,2)
    %         if sum(AUC(1,i)>AUCit(:,i))>fix(sig_threshold*size(AUCit,1))
    %             psig_AUC(i,subs)=1;
    %         else
    %             psig_AUC(i,subs)=0;
    %         end
    %     end
    
    %% Preparing for Plotting
    clearvars X Xp
    cond1=1;
    cond2=2;
    cond3=3;
    cond4=4;
    c=0;
    for t=tws
        c=c+1;
        X(:,:,c)=[mean(signals{1,cond1}(:,t,:),3) mean(signals{1,cond2}(:,t,:),3) mean(signals{1,cond3}(:,t,:),3) mean(signals{1,cond4}(:,t,:),3);1 1 1 1];
        Xp(:,c)=squeeze(wp(:,indx_max_AUC(subs))'*squeeze(X(:,:,c)));
    end
    XXX(:,:,:,subs)=X;
    XXXp(:,:,subs)=Xp;
    
    % Scalp
    a(:,subs)=(squeeze(X(:,:,indx_max_AUC(subs)))*squeeze(Xp(:,indx_max_AUC(subs))))./(squeeze(Xp(:,indx_max_AUC(subs)))'*squeeze(Xp(:,indx_max_AUC(subs))));
    
    % Matrix format sorted trials    
    Xp_sorted=nan*ones(240,(2000-span)./steps);
    if cats==0
        cond_size=[60 60 60 60]*2;
    else
        cond_size=[120 40 40 40]*2;        
    end
        
    for cond=1:4
        clearvars XX Xpp
        c=0;
        for t=tws
            c=c+1;
            XX=[squeeze(signals{1,cond}(:,t,:));ones(1,size(signals{1,cond},3))];
            Xpp(:,c)=squeeze(wp(:,indx_max_AUC(subs))'*XX);
        end
        [~,indx_Xp_sorted]=sort(Xpp(:,indx_max_AUC(subs)),'ascend');
        Xp_sorted(sum(cond_size(1:cond))-length(indx_Xp_sorted)+1:sum(cond_size(1:cond)),:)=Xpp(indx_Xp_sorted,:);
    end
    Xp_sorted_t(:,:,subs)=Xp_sorted; 
    clearvars wp wt wtb    
    subs
    
end
%% Plotting
subjects=[1:length(subjects)];

if stim_resp==1
    span=[-500:steps:1449];
else
    span=[-1500:steps:449];
end

figure;
plot(span,nanmean(AUCn(subjects,:),1))
hold on;
plot(span,nanmean(AUC(subjects,:),1))
plot(span,nanmean(AUC_after(subjects,:),1))
legend ('Before transformation','After trans W across time','After trans W best time')
xlabel('Time [ms]')
ylabel('AUC')
line([span(1) span(end)],[0.5 0.5])
line([0 0],[0.4 1])
if stim_resp==1
    xlim([-500 1000])
else
    xlim([-1000 500])
end

figure;
plot(span,squeeze(nanmean(nanmean(XXX(:,1,:,subjects)),4))); hold on;
plot(span,squeeze(nanmean(nanmean(XXX(:,2,:,subjects)),4)));
plot(span,squeeze(nanmean(nanmean(XXX(:,3,:,subjects)),4)));
plot(span,squeeze(nanmean(nanmean(XXX(:,4,:,subjects)),4)));
% plot(abs(squeeze(mean(mean(XXX(:,1,:,:)),4))-squeeze(mean(mean(XXX(:,4,:,:)),4))));
if cats==1
    legend ('Control','Famous','Familiar','Self','Location','northwest');
else
    legend ('Coherence = 0.22','Coherence = 0.33','Coherence = 0.45','Coherence = 0.55','Location','northwest');
end
xlabel('Time [ms]')
ylabel('Amplitude [uv]')
hold on;
line([0 0],[-4 6])
if stim_resp==1
    xlim([-500 1000])
else
    xlim([-1000 500])
end

figure;
plot(span,squeeze(nanmean(XXXp(:,:,subjects),3))'); hold on;
% plot(abs(squeeze(mean(XXXp(1,:,:),3))-squeeze(mean(XXXp(4,:,:),3))));
if cats==1
    legend ('Control','Famous','Familiar','Self','Location','northwest');
else
    legend ('Coherence = 0.22','Coherence = 0.33','Coherence = 0.45','Coherence = 0.55','Location','northwest');
end
xlabel('Time [ms]')
ylabel('Mean y [uv]')
hold on;
line([0 0],[0 0.2])
if stim_resp==1
    xlim([-500 1000])
else
    xlim([-1000 500])
end


figure;
% Xp_sorted_t(isoutlier(Xp_sorted_t,1))=nan;
tmp=nanmean(Xp_sorted_t(:,:,subjects),3);
imagesc(tmp(end:-1:1,:))
if stim_resp==1
    imagesc(tmp(end:-1:1,1:end-45))
else
    imagesc(tmp(end:-1:1,50:end))
end

colorbar

figure;
if stim_resp==1
    histogram(indx_max_AUC.*steps-500,20)
    xlim([-500 1000])
else
    histogram(indx_max_AUC.*steps-1500,20)
    xlim([-1000 500])
end