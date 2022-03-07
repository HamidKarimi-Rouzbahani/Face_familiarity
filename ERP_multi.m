clc;
close all;
clear all;
load('All_subj_preprocessed_datasets.mat');

%% Settings
cats=0; % categories=1 and noise_levels=0
only_answered_trials=1; % all trials=0
only_correct_trials=1; % plus incorrect=0
stim_resp=2; % for stim=1, for resp=2

if stim_resp==1
    epoch_span=[-200 1499]; %for stimulus aligned
else
    epoch_span=[-1499 200]; %for response aligned
end

subjs=[1:16 20 21];      % subj #

%% Behavioral
RT_cat=nan*ones(20,4);
RT_noise=nan*ones(20,4);
Acc_cat=nan*ones(20,4);
Acc_noise=nan*ones(20,4);
for subj=subjs
    nanmean(Exp_data{1, subj}.stim.ResponseData.Values(2,:))
    for i=1:4
        ids=(Exp_data{1,subj}.stim.ResponseData.Values(5,:)==i);
        if only_correct_trials==1
            ids=ids(Exp_data{1,subj}.stim.ResponseData.Values(2,:)==1);
        end
          RT_cat(subj,i)=round(nanmean(Exp_data{1,subj}.stim.ResponseData.Values(1,ids)));
          Acc_cat(subj,i)=nanmean(Exp_data{1,subj}.stim.ResponseData.Values(2,ids));

        noises=[0.22 0.3 0.45 0.55];
        ids=(Exp_data{1,subj}.stim.ResponseData.Values(4,:)==noises(i));
        if only_correct_trials==1
            ids=ids(Exp_data{1,subj}.stim.ResponseData.Values(2,:)==1);
        end
          RT_noise(subj,i)=round(nanmean(Exp_data{1,subj}.stim.ResponseData.Values(1,ids)));  
          Acc_noise(subj,i)=nanmean(Exp_data{1,subj}.stim.ResponseData.Values(2,ids));  
    end
%     xlabel('Reaction time [ms]');
%     ylabel('Number of trials [n]');
%     legend data MeanOfAll Control Famous Familiar Self 0.22 0.3 0.45 0.55
end
errorbar(nanmean(RT_cat),nanstd(RT_cat)./sqrt(length(subjs)))
hold on;
errorbar(nanmean(RT_noise),nanstd(RT_cat)./sqrt(length(subjs)))
legend Categories Coherence
figure;
errorbar(nanmean(Acc_cat),nanstd(Acc_cat)./sqrt(length(subjs)))
hold on;
errorbar(nanmean(Acc_noise),nanstd(Acc_noise)./sqrt(length(subjs)))
legend Categories Coherence
ccc

%% ERPs
for ch=0:8:56
    % for ch=0:1
    figure;
    c=0;
    for chf=1+ch:8+ch
        c=c+1;
        cc=0;
        %         indx=1:size(EEG_signals{stim_resp,subj},3); %all
        if cats==1
            for counte=[1 2 4 3]
                colors=[0 0.45 0.7 1];
                cc=cc+1;
                signals=nan.*ones(20,1700);
                for subj=subjs
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
                    signals(subj,:)=squeeze(mean(EEG_signals{stim_resp,subj}(chf,:,indx{cc,subj}),3));
                end
                subplot(2,4,c);
                plot([epoch_span(1):epoch_span(2)],nanmean(signals),'Color',[0.8 0.2 1]*colors(cc),'LineWidth',1.5);
                hold on;
            end
        else
            for counte=[0.22 0.3 0.45 0.55]
                colors=[0 0.45 0.7 1];
                cc=cc+1;
                signals=nan.*ones(20,1700);
                for subj=subjs
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
                    signals(subj,:)=squeeze(mean(EEG_signals{stim_resp,subj}(chf,:,indx{cc,subj}),3));
                end
                subplot(2,4,c);
                plot([epoch_span(1):epoch_span(2)],nanmean(signals),'Color',[1 0.2 0.8]*colors(cc),'LineWidth',1.5);
                hold on;
            end
        end
        line([epoch_span(1) epoch_span(2)],[0 0]);
        line([0 0],[-10 10]);
        xlim([epoch_span(1) epoch_span(2)]);
        %         ylim([-20 20]);
        xlabel('Time [ms]');
        ylabel('Amplitude [uv]');
        title(chan_locations{1,1}{chf,1})
        hold off;
    end
    if cats==1
        legend ('Control','Famous','Familiar','Self','Location','southeast');
    else
        legend ('Noise = 0.22','Noise = 0.30','Noise = 0.45','Noise = 0.55','Location','southeast');
    end
end
cccc
%% Whole brain ERPs
chf=[1:64];
cc=0;
subjs=[1 3 5];      % subj #
for counte1=[0.22 0.3 0.45 0.55]
    cc=cc+1;
    dd=0;
    for counte2=[1 2 4 3]
        dd=dd+1;
        signals=nan.*ones(20,1700);
        for subj=subjs
            indx{cc,dd,subj}=(([Exp_data{1,subj}.stim.stimTrain.imageNoise]==counte1) & ([Exp_data{1,subj}.stim.stimTrain.imageCategory]==counte2));
            if only_answered_trials==1
                tmp=indx{cc,dd,subj};
                tmp(isnan(Exp_data{1,subj}.stim.ResponseData.Values(2,:)))=0;
                indx{cc,dd,subj}=tmp;
                if only_correct_trials==1
                    tmp=indx{cc,dd,subj};
                    tmp(Exp_data{1,subj}.stim.ResponseData.Values(2,:)==0)=0;
                    indx{cc,dd,subj}=tmp;
                end
            end
            signals(subj,:)=squeeze(nanmean(nanmean(EEG_signals{stim_resp,subj}(chf,:,indx{cc,dd,subj}),3),1));
        end
        subplot(2,2,cc);
        plot([epoch_span(1):epoch_span(2)],nanmean(signals),'LineWidth',1.5);
        title(num2str(counte1));
        hold on;
    end
    line([epoch_span(1) epoch_span(2)],[0 0]);
    line([0 0],[-10 10]);
    xlim([epoch_span(1) epoch_span(2)]);
    %     ylim([-20 20]);
    xlabel('Time [ms]');
    ylabel('Amplitude [uv]');
    hold off;
end
legend ('Control','Famous','Familiar','Self','Location','southeast');

