clc;
close all;
clear all;
%% Behavioral
% subjs=[1:16 20 21];    % subj #
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','20','21'};
TPR=nan*ones(21,5,4,2);
FPR=nan*ones(21,5,4,2);
RT_correct=nan*ones(21,5,4,2);
c=0;
coherences=[0.22 0.3 0.45 0.55];
% for subj=13:length(subject)
for subj=1:length(subject)
    c=c+1;
    load(['Subject_',subject{c},'_preprosessed.mat'],'Exp_data','chan_locations');
    
    for coherence=1:4
        if subj==13
            sessions=1;
        else
            sessions=[1 2];
        end
        for session=sessions
            if isempty(Exp_data{1,session})==0
                
                ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]==1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_1=tmp;
                
                ids_2=([Exp_data{1,session}.stim.stimTrain.imageCategory]>1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                tmp=ids_2;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_2=tmp;
                               
                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==0 | Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_3=tmp;

                TPR(subj,1,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
                FPR(subj,1,coherence,session)=1-(nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_2)));
                RT_correct(subj,1,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(1,ids_3));

                             
                ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]>1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_1=tmp;
                
                ids_2=([Exp_data{1,session}.stim.stimTrain.imageCategory]==1 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                tmp=ids_2;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_2=tmp;
                               
                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==0 | Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_3=tmp;

                TPR(subj,2,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
                FPR(subj,2,coherence,session)=1-(nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_2)));
                RT_correct(subj,2,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(1,ids_3));

                % subfamiliarity levels
                ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]==2 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_1=tmp;

                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==0 | Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_3=tmp;

                TPR(subj,3,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
                RT_correct(subj,3,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(1,ids_3));

                % subfamiliarity levels
                ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]==4 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_1=tmp;

                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==0 | Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_3=tmp;

                TPR(subj,4,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
                RT_correct(subj,4,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(1,ids_3));

                
                ids_1=([Exp_data{1,session}.stim.stimTrain.imageCategory]==3 & [Exp_data{1, session}.stim.stimTrain.imageNoise]==coherences(coherence));
                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_1=tmp;

                tmp=ids_1;
                tmp([Exp_data{1,session}.stim.ResponseData.Values(2,:)==0 | Exp_data{1,session}.stim.ResponseData.Values(2,:)==inf])=0;
                ids_3=tmp;

                TPR(subj,5,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(2,ids_1));
                RT_correct(subj,5,coherence,session)=nanmean(Exp_data{1,session}.stim.ResponseData.Values(1,ids_3));
                
                
            end
        end
    end
end

plot(squeeze(nanmean(nanmean(TPR,1),4))','LineWidth',3);
xlabel('Coherence')
ylabel ('TPR')
legend Control Familiar Famous Personally Self

figure;
plot(squeeze(nanmean(nanmean(FPR,1),4))','LineWidth',3);
xlabel('Coherence')
ylabel ('FPR')
legend Control Familiar

figure;
plot(squeeze(nanmean(nanmean(RT_correct,1),4))','LineWidth',3);
xlabel('Coherence')
ylabel ('RT')
legend Control Familiar Famous Personally Self

save('Hit_rates_False_alarms_and_Correct_reaction_times.mat','TPR','FPR','RT_correct')
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
%% Dprime and bias
hit=squeeze(nanmean(nanmean(nanmean(TPR(:,1,:,1),1),2),4));
falseA=squeeze(nanmean(nanmean(nanmean(FPR(:,1,:,1),1),2),4));

[dpri,ccrit] = dprime(hit,falseA)

