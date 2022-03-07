clc;
close all;
clear all;
%% Behavioral
% subjs=[1:16 20 21];    % subj #
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','20','21'};
RT_cat=nan*ones(21,4,4,2);
RT_noise=nan*ones(21,4,4,2);
Acc_cat=nan*ones(21,4,4,2);
Acc_noise=nan*ones(21,4,4,2);
c=0;
coherences=[0.22 0.3 0.45 0.55];
for subj=1:length(subject)
    c=c+1;
    load(['Subject_',subject{c},'_preprosessed.mat'],'Exp_data','chan_locations');
    for coherence=1:4
        if subj~=3
            sessions=[1 2];
        else
            sessions=1;
        end
        for session=sessions
            if isempty(Exp_data{1,session})==0
                %                     ids=(Exp_data{1,session}.stim.ResponseData.Values(5,:)==cat);
                
                ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]==1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                ids_2=([Exp_data{1,session}.stim.stimTrain.imageCategory]>1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                TPR(subj,1,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
                FPR(subj,1,coherence,session)=1-(nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_2)));


                ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]==2 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                ids_2=([Exp_data{1,session}.stim.stimTrain.imageCategory]==1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                TPR(subj,2,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
                FPR(subj,2,coherence,session)=(1-(nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_2))));

                ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]==3 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                ids_2=([Exp_data{1,session}.stim.stimTrain.imageCategory]==1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                TPR(subj,3,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
                FPR(subj,3,coherence,session)=(1-(nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_2))));
                
                ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]==4 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                ids_2=([Exp_data{1,session}.stim.stimTrain.imageCategory]==1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                TPR(subj,4,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
                FPR(subj,4,coherence,session)=(1-(nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_2))));

%                 ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]>1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
%                 ids_2=([Exp_data{1,session}.stim.stimTrain.imageCategory]==1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
%                 TPR(subj,5,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
%                 FPR(subj,5,coherence,session)=1-(nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_2)));

                for cat=1:4
                    dprimes(subj,cat,coherence,session)=dprime(TPR(subj,cat,coherence,session),FPR(subj,cat,coherence,session));
                end
                %                     Acc_cat(subj,cat,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids));
                %                     RT_cat(subj,cat,coherence,session)=round(nanmean(Exp_data{1,session}.stim.ResponseData.Values(1,ids)));
            end
        end
    end
end

plot(squeeze(nanmean(nanmean(dprimes,1),4))','LineWidth',3);
xlabel('Coherence')
ylabel ('dprime')
% xticklable({'22','33','45','55'})
legend Control Famous Self Personallyfamiliar

% gca()

% errorbar(nanmean(nanmean(RT_cat),3),nanstd(nanmean(RT_cat,3))./sqrt(length(subject)))
% hold on;
% errorbar(nanmean(nanmean(RT_noise),3),nanstd(nanmean(RT_noise,3))./sqrt(length(subject)))
% xlabel('Conditions');
% ylabel('Reaction time [ms]');
% legend Categories Coherence
% figure;
% errorbar(nanmean(nanmean(Acc_cat),3),nanstd(nanmean(Acc_cat,3))./sqrt(length(subject)))
% hold on;
% errorbar(nanmean(nanmean(Acc_noise),3),nanstd(nanmean(Acc_noise,3))./sqrt(length(subject)))
% xlabel('Conditions');
% ylabel('Proportion correct');
% legend Categories Coherence
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
