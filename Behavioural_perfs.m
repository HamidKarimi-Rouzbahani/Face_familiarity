clc;
close all;
clear all;
% load('./Session1/All_subj_preprocessed_datasets_baseline_intact_session1.mat');
% subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','20','21'};
% c=0;
% for i=[1:16 20 21]
%     c=c+1
%     EEG_signal{1,1}=EEG_signals{1,i};
%     EEG_signal{2,1}=EEG_signals{2,i};
%     Exp_datum{1,1}=Exp_data{1,i};
%     save(['Subject_',subject{c},'_preprosessed_signals_session1.mat'],'EEG_signal','Exp_datum','chan_locations');
% end
% cc
% load('./Session2/All_subj_preprocessed_datasets_baseline_intact_session2.mat');
% c=0;
% for i=[1:16 20 21]
%     c=c+1
%     EEG_signal{1,1}=EEG_signals{1,i};
%     EEG_signal{2,1}=EEG_signals{2,i};
%     Exp_datum{1,1}=Exp_data{1,i};
%     save(['Subject_',subject{c},'_preprosessed_signals_session2.mat'],'EEG_signal','Exp_datum','chan_locations');
% end
%
% d=0;
% for i=[1:16 20 21]
%     d=d+1
%     load(['Subject_',subject{d},'_preprosessed_signals_session1.mat']);
%     EEG_signals{1,1}=EEG_signal{1,1};
%     EEG_signals{2,1}=EEG_signal{2,1};
%     Exp_data{1,1}=Exp_datum{1,1};
%     load(['Subject_',subject{d},'_preprosessed_signals_session2.mat']);
%     EEG_signals{1,2}=EEG_signal{1,1};
%     EEG_signals{2,2}=EEG_signal{2,1};
%     Exp_data{1,2}=Exp_datum{1,1};
%     save(['Subject_',subject{d},'_preprosessed.mat'],'EEG_signals','Exp_data','chan_locations','-v7.3');
% end

%% Behavioral
only_correct_trials=0; % plus incorrect=0
% subjs=[1:16 20 21];    % subj #
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','20','21'};
RT_cat=nan*ones(21,4,2);
RT_noise=nan*ones(21,4,2);
Acc_cat=nan*ones(21,4,2);
Acc_noise=nan*ones(21,4,2);
c=0;
for subj=1:length(subject)
    c=c+1;
    load(['Subject_',subject{c},'_preprosessed.mat'],'Exp_data','chan_locations');
    for cat=1:4
        for cat=1:4
            for session=1:2
                if isempty(Exp_data{1,session})==0
                    ids=(Exp_data{1,session}.stim.ResponseData.Values(5,:)==cat);
                    if only_correct_trials==1
                        ids=ids(Exp_data{1,session}.stim.ResponseData.Values(2,:)==1);
                    end
                    RT_cat(subj,cat,session)=round(nanmean(Exp_data{1,session}.stim.ResponseData.Values(1,ids)));
                    Acc_cat(subj,cat,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids));
                    
                    %                 noises=[0.22 0.3 0.45 0.55];
                    %                 ids=(Exp_data{1,session}.stim.ResponseData.Values(4,:)==noises(i));
                    %                 Acc_noise(subj,i,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids));
                    %                 if only_correct_trials==1
                    %                     ids=ids(Exp_data{1,session}.stim.ResponseData.Values(2,:)==1);
                    %                 end
                    %                 RT_noise(subj,i,session)=round(nanmean(Exp_data{1,session}.stim.ResponseData.Values(1,ids)));
                end
            end
        end
    end
    %     legend data MeanOfAll Control Famous Familiar Self 0.22 0.3 0.45 0.55
end

errorbar(nanmean(nanmean(RT_cat),3),nanstd(nanmean(RT_cat,3))./sqrt(length(subject)))
hold on;
errorbar(nanmean(nanmean(RT_noise),3),nanstd(nanmean(RT_noise,3))./sqrt(length(subject)))
xlabel('Conditions');
ylabel('Reaction time [ms]');
legend Categories Coherence
figure;
errorbar(nanmean(nanmean(Acc_cat),3),nanstd(nanmean(Acc_cat,3))./sqrt(length(subject)))
hold on;
errorbar(nanmean(nanmean(Acc_noise),3),nanstd(nanmean(Acc_noise,3))./sqrt(length(subject)))
xlabel('Conditions');
ylabel('Proportion correct');
legend Categories Coherence
ccc
%% Significance
for cat=1:4
    for j=1:4
        significance_RT_cat(cat,j)=signrank(nanmean(RT_cat(:,cat,:),3),nanmean(RT_cat(:,j,:),3));
    end
end

for cat=1:4
    for j=1:4
        significance_RT_noise(cat,j)=signrank(nanmean(RT_noise(:,cat,:),3),nanmean(RT_noise(:,j,:),3));
    end
end

for cat=1:4
    for j=1:4
        significance_Acc_cat(cat,j)=signrank(nanmean(Acc_cat(:,cat,:),3),nanmean(Acc_cat(:,j,:),3));
    end
end

for cat=1:4
    for j=1:4
        significance_Acc_noise(cat,j)=signrank(nanmean(Acc_noise(:,cat,:),3),nanmean(Acc_noise(:,j,:),3));
    end
end
