clc;
clear all;
% close all;
%% Loading and analysis of weights
only_answered_trials=1; % all trials=0
only_correct_trials=1; % plus incorrect=0
cats=1;
stim_resp=2;
subs=0;
subjects=[1:16 20 21];
Xp_sorted_t=nan*ones(480,1950,length(subjects));
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
sessions=[1 2];
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
    %% EOG AND M1 are channels #32 and #13
    for i=1:4
        signals{1,i}(13,:,:)=[];
        signals{1,i}(32,:,:)=[];
    end
    clearvars -except subject sessions Exp_data indx_max_AUC Xp_sorted XXX XXXp AUC AUCn a subs EEG_signals signals subj cats stim_resp only_answered_trials only_correct_trials
    %%
    load(['weights_subj_',num2str(subj),'diff_tm.mat']);
    cond1=1;
    cond2=4;
    span=50;
    tws=1:2000-span;
    clearvars X
    c=0;
    d=[ones(1,size(signals{1,cond1}(:,1,:),3)) zeros(1,size(signals{1,cond2}(:,1,:),3))];
    for t=tws
        c=c+1;
        X=[mean(signals{1,cond1}(:,t,:),3) mean(signals{1,cond2}(:,t,:),3); 1 1];
        Xp=wp'*squeeze(X);
        Xpn=mean(X);
        XX=[squeeze(signals{1,cond1}(:,t,:)) squeeze(signals{1,cond2}(:,t,:));ones(1,size(signals{1,cond1}(:,t,:),3)+size(signals{1,cond2}(:,t,:),3))];
        XXp=wp'*XX;
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
    [~,indx_max_AUC(subs)]=max(AUC(subs,:));    
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
        Xp(:,c)=squeeze(wp'*squeeze(X(:,:,c)));
    end
    XXX(:,:,:,subs)=X;
    XXXp(:,:,subs)=Xp;
    
    % Scalp
    if stim_resp==1
        training_time=450;
    else
        training_time=1370;
    end
    a(:,subs)=(squeeze(X(:,:,training_time))*squeeze(Xp(:,training_time)))./(squeeze(Xp(:,training_time))'*squeeze(Xp(:,training_time)));
    
    % Matrix format sorted trials
    Xp_sorted=nan*ones(480,2000-span);
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
            Xpp(:,c)=wp'*XX;
        end
        [~,indx_Xp_sorted]=sort(Xpp(:,training_time),'ascend');
        Xp_sorted(sum(cond_size(1:cond))-length(indx_Xp_sorted)+1:sum(cond_size(1:cond)),:)=Xpp(indx_Xp_sorted,:);
    end
    Xp_sorted_t(:,:,subs)=Xp_sorted;
    clearvars wp wt wtb
    subs    
end

%% Plotting
subjects=[1:18];


if stim_resp==1
    span=[-500:1449];
else
    span=[-1500:449];
end

figure;
plot(span,nanmean(AUCn(subjects,:),1))
hold on;
plot(span,nanmean(AUC(subjects,:),1))
legend ('Before transformation','After transformation')
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
    legend ('Coherence = 0.22','Coherence = 0.30','Coherence = 0.45','Coherence = 0.55','Location','northwest');
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
Xp_sorted_t(isoutlier(Xp_sorted_t,1))=nan;
tmp=nanmean(Xp_sorted_t(:,:,subjects),3);
imagesc(tmp(end:-1:1,:))
colorbar

