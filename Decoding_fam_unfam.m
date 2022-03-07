clc;
clear all;
close all;

%% Settings
% cats=0; % categories=1 and coh=0
only_answered_trials=1; % all trials=0
only_correct_trials=1; % plus incorrect=0 you should also change the saving file names
only_incorrect_trials=0; % plus correct=0

only_answered_trials2=1; % all trials=0
only_correct_trials2=0; % plus incorrect=0 you should also change the saving file names
only_incorrect_trials2=1; % plus correct=0
% stim_resp=1; % for stim=1, for resp=2
baseline_correction=1;
steps=10;

%% ERPs
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
sessions=[1 2];
accuracy=nan*ones(21,6,200);
test_on_errors=0;
for stim_resp=[1 2]
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
        
        % EOG AND M1 are channels #32 and #13
        for i=1:4
            signals{1,i}(13,:,:)=[];
            signals{1,i}(32,:,:)=[];
        end
        % Adding all familiar categories to cond #2
        sizes_signals=[size(signals{1,1},3) size(signals{1,2},3) size(signals{1,3},3) size(signals{1,4},3)];
        signals{1,2}(:,:,end+1:end+size(signals{1,3},3))=signals{1,3};
        signals{1,2}(:,:,end+1:end+size(signals{1,4},3))=signals{1,4};
        signals{1,3}=[];
        signals{1,4}=[];
        clearvars -except test_on_errors sizes_signals accuracy Exp_data baseline_correction EEG_signals subject sessions steps baseline_span signals2 signals subj cats stim_resp only_answered_trials only_correct_trials only_incorrect_trials only_answered_trials2 only_correct_trials2 only_incorrect_trials2
        
        if test_on_errors==1
            %% Second dataset
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
                    indx{cc,session}=[Exp_data{1,session}.stim.stimTrain.imageCategory]==counte;
                    if only_answered_trials2==1
                        tmp=indx{cc,session};
                        tmp(isnan(Exp_data{1,session}.stim.ResponseData.Values(2,:)))=0;
                        indx{cc,session}=tmp;
                        if only_correct_trials2==1
                            tmp=indx{cc,session};
                            tmp(Exp_data{1,session}.stim.ResponseData.Values(2,:)==0)=0;
                            indx{cc,session}=tmp;
                        end
                        if only_incorrect_trials2==1
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
                signals2{cc}=tmpt;
            end
            
            
            %             if baseline_correction
            %                 for i=1:4
            %                     signals2{1,i}=signals2{1,i}-repmat(nanmean(nanmean(signals2{1,i}(:,baseline_span,:),3),2),[1 2000 size(signals2{1,i},3)]);
            %                 end
            %             end
            
            if baseline_correction
                tt=nanmean([nanmean(nanmean(signals2{1,1}(:,baseline_span,:),3),2) nanmean(nanmean(signals2{1,2}(:,baseline_span,:),3),2) nanmean(nanmean(signals2{1,3}(:,baseline_span,:),3),2) nanmean(nanmean(signals2{1,4}(:,baseline_span,:),3),2)],2);
                for i=1:4
                    signals2{1,i}=signals2{1,i}-repmat(tt,[1 2000 size(signals2{1,i},3)]);
                end
            end
            
            % if cats==1
            %     legend ('Control','Famous','Familiar','Self','Location','southeast');
            % else
            %     legend ('Coherence = 0.22','Coherence = 0.30','Coherence = 0.45','Coherence = 0.55','Location','southeast');
            % end
            
            % EOG AND M1 are channels #32 and #13
            for i=1:4
                signals2{1,i}(13,:,:)=[];
                signals2{1,i}(32,:,:)=[];
            end
            
            % Adding all familiar categories to cond #2
            sizes_signals=[size(signals2{1,1},3) size(signals2{1,2},3) size(signals2{1,3},3) size(signals2{1,4},3)];
            signals2{1,2}(:,:,end+1:end+size(signals2{1,3},3))=signals2{1,3};
            signals2{1,2}(:,:,end+1:end+size(signals2{1,4},3))=signals2{1,4};
            signals2{1,3}=[];
            signals2{1,4}=[];
            clearvars -except test_on_errors sizes_signals accuracy subject sessions steps baseline_span baseline_correction signals signals2 subj cats stim_resp only_answered_trials only_correct_trials only_incorrect_trials only_answered_trials2 only_correct_trials2 only_incorrect_trials2
            
            %% SVM
            cond1=1;
            cond2=2;
            c=0;
            for tws=1:steps:2000-steps
                tw=tws:tws+steps;
                c=c+1;
                X1=[squeeze(mean(signals{1,cond1}(:,tw,:),2))';squeeze(mean(signals{1,cond2}(:,tw,:),2))'];
                X2=[squeeze(mean(signals2{1,cond1}(:,tw,:),2))';squeeze(mean(signals2{1,cond2}(:,tw,:),2))'];
                y1=[ones(size(squeeze(signals{1,cond1}(:,tws,:))',1),1);zeros(size(squeeze(signals{1,cond2}(:,tws,:))',1),1)];
                y2=[ones(size(squeeze(signals2{1,cond1}(:,tws,:))',1),1);zeros(size(squeeze(signals2{1,cond2}(:,tws,:))',1),1)];
                
                SVMModel = fitcsvm(X1,y1);
                
                [label,~] = predict(SVMModel,X2);
                
                % Unfamiliar vs Familiar
                accuracy(subj,1,c)=sum(label==y2)./length(label);
                
                % Unfamiliar and Familiar separately
                accuracy(subj,2,c)=sum(label(1:size(squeeze(signals2{1,cond1}(:,tws,:))',1))==y2(1:size(squeeze(signals2{1,cond1}(:,tws,:))',1)))./size(squeeze(signals2{1,cond1}(:,tws,:))',1);
                accuracy(subj,3,c)=sum(label(size(squeeze(signals2{1,cond1}(:,tws,:))',1)+1:end)==y2(size(squeeze(signals2{1,cond1}(:,tws,:))',1)+1:end))./size(squeeze(signals2{1,cond2}(:,tws,:))',1);
                
                % Familiar categories in order: Famous, Familiar, Self
                accuracy(subj,4,c)=sum(label(sizes_signals(1)+1:sizes_signals(1)+sizes_signals(2))==y2(sizes_signals(1)+1:sizes_signals(1)+sizes_signals(2)))./sizes_signals(2);
                accuracy(subj,5,c)=sum(label(sizes_signals(1)+sizes_signals(2)+1:sizes_signals(1)+sizes_signals(2)+sizes_signals(3))==y2(sizes_signals(1)+sizes_signals(2)+1:sizes_signals(1)+sizes_signals(2)+sizes_signals(3)))./sizes_signals(3);
                accuracy(subj,6,c)=sum(label(sizes_signals(1)+sizes_signals(2)+sizes_signals(3)+1:sizes_signals(1)+sizes_signals(2)+sizes_signals(3)+sizes_signals(4))==y2(sizes_signals(1)+sizes_signals(2)+sizes_signals(3)+1:sizes_signals(1)+sizes_signals(2)+sizes_signals(3)+sizes_signals(4)))./sizes_signals(4);
                
                [test_on_errors stim_resp subj c]
            end
            if stim_resp==1
                save(['st_Decoding_fam_unfam_tr_cor_ts_incor_baselined_windowed2.mat'],'accuracy');
            else
                save(['rp_Decoding_fam_unfam_tr_cor_ts_incor_baselined_windowed2.mat'],'accuracy');
            end
        else
            %% SVM
            cond1=1;
            cond2=2;
            c=0;
            for tws=1:steps:2000-steps
                tw=tws:tws+steps;
                c=c+1;
                X=[squeeze(mean(signals{1,cond1}(:,tw,:),2))';squeeze(mean(signals{1,cond2}(:,tw,:),2))'];
                y=[ones(size(squeeze(signals{1,cond1}(:,tws,:))',1),1);zeros(size(squeeze(signals{1,cond2}(:,tws,:))',1),1)];
                for iter=1:10
                    tr_inds=randsample([1:size(X,1)],fix(size(X,1)*0.8));
                    ts_inds=[1:size(X,1)];
                    ts_inds(tr_inds)=[];
                    X1=X(tr_inds,:);
                    y1=y(tr_inds);
                    SVMModel = fitcsvm(X1,y1);
                    X2=X(ts_inds,:);
                    y2=y(ts_inds);
                    [label,~] = predict(SVMModel,X2);
                    
                    inds_tst1=ones(1,length(ts_inds));
                    inds_tst2=(ts_inds<sizes_signals(1)+1);
                    inds_tst3=(ts_inds>sizes_signals(1));
                    inds_tst4=(ts_inds>sizes_signals(1) & ts_inds<sizes_signals(1)+sizes_signals(2)+1);
                    inds_tst5=(ts_inds>sizes_signals(1)+sizes_signals(2) & ts_inds<sizes_signals(1)+sizes_signals(2)+sizes_signals(3)+1);
                    inds_tst6=(ts_inds>sizes_signals(1)+sizes_signals(2)+sizes_signals(3)+1 & ts_inds<sizes_signals(1)+sizes_signals(2)+sizes_signals(3)+sizes_signals(4)+1);
                    
                    % Unfamiliar vs Familiar
                    tmp(1,iter)=sum(label(inds_tst1)==y2(inds_tst1))./sum(inds_tst1);
                    
                    % Unfamiliar and Familiar separately
                    tmp(2,iter)=sum(label(inds_tst2)==y2(inds_tst2))./sum(inds_tst2);
                    tmp(3,iter)=sum(label(inds_tst3)==y2(inds_tst3))./sum(inds_tst3);
                    
                    % Familiar categories in order: Famous, Familiar, Self
                    tmp(4,iter)=sum(label(inds_tst4)==y2(inds_tst4))./sum(inds_tst4);
                    tmp(5,iter)=sum(label(inds_tst5)==y2(inds_tst5))./sum(inds_tst5);
                    tmp(6,iter)=sum(label(inds_tst6)==y2(inds_tst6))./sum(inds_tst6);
                    [test_on_errors stim_resp subj c iter]
                end
                accuracy(subj,:,c)=squeeze(nanmean(tmp,2));
            end
            if stim_resp==1
                save(['st_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed2.mat'],'accuracy');
            else
                save(['rp_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed2.mat'],'accuracy');
            end
        end
    end
end
%% plotting
stim_resp=1;
baseline_correction=1;
windowed=1;
cor_cor=1;


if baseline_correction==0
    if cor_cor==1
        if stim_resp==1
            load(['st_Decoding_fam_unfam_tr_cor_ts_cor.mat']);
        else
            load(['rp_Decoding_fam_unfam_tr_cor_ts_cor.mat']);
        end
    else
        if stim_resp==1
            load('st_Decoding_fam_unfam_tr_cor_ts_incor.mat');
        else
            load('rp_Decoding_fam_unfam_tr_cor_ts_incor.mat');
        end
    end
else
    if windowed==0
        if cor_cor==1
            if stim_resp==1
                load(['st_Decoding_fam_unfam_tr_cor_ts_cor_baselined.mat']);
            else
                load(['rp_Decoding_fam_unfam_tr_cor_ts_cor_baselined.mat']);
            end
        else
            if stim_resp==1
                load('st_Decoding_fam_unfam_tr_cor_ts_incor_baselined.mat');
            else
                load('rp_Decoding_fam_unfam_tr_cor_ts_incor_baselined.mat');
            end
        end
    else
        if cor_cor==1
            if stim_resp==1
                load(['st_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed2.mat']); % windowed= incorrectly baselined, windowed2 is correct
            else
                load(['rp_Decoding_fam_unfam_tr_cor_ts_cor_baselined_windowed2.mat']);
            end
        else
            if stim_resp==1
                load('st_Decoding_fam_unfam_tr_cor_ts_incor_baselined_windowed2.mat');
            else
                load('rp_Decoding_fam_unfam_tr_cor_ts_incor_baselined_windowed2.mat');
            end
        end
    end
end
steps=10;
if stim_resp==1
    span=[-500:steps:1490];
else
    span=[-1500:steps:490];
end
figure;
for i=1:6
    plot(span,smooth(squeeze(nanmean(accuracy([1:9 12:21],i,:),1)),10),'linewidth',3);
    hold on;
end
line([span(1) span(end)],[0.5 0.5])
line([0 0],[0.2 0.7])
legend ('Familiar vs Unfamiliar','Unfamiliar','Familiar','Famous','Personally familiar','Self','Location','northwest');
xlabel('Time [ms]')
ylabel('Decoding accuracy')
hold on;
% if stim_resp==1
%     xlim([-500 1000])
% else
%     xlim([-1000 500])
% end

