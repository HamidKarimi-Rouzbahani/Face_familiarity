% clc;
% close all;
% clear all;

%% Settings
Subj='03';
load(['Face_Discrimination_Data_Subject_',Subj,'_0000.mat']);epoch_span=[-200 1028]; %for stimulus aligned
% epoch_span=[-1028 200]; %for response aligned
cats=1; % categories=1 and noise_levels=0
only_answered_trials=1; % all trials=0
only_correct_trials=1; % plus incorrect=0
ds=6;

%% Behavioral
histogram(stim.ResponseData.Values(1,:),40)
hold on;
line([round(nanmean(stim.ResponseData.Values(1,:))) round(nanmean(stim.ResponseData.Values(1,:)))],[0 20],'LineWidth',2)
for i=1:4
    ids=(stim.ResponseData.Values(5,:)==i);
    if only_correct_trials==1
        ids=ids(stim.ResponseData.Values(2,:)==1);
    end
    line([round(nanmean(stim.ResponseData.Values(1,ids))) round(nanmean(stim.ResponseData.Values(1,ids)))],[0 20],'Color',[0.8 0 0]*i./4,'LineWidth',2);
    
    noises=[0.22 0.3 0.45 0.55];
    ids=(stim.ResponseData.Values(4,:)==noises(i));
    if only_correct_trials==1
        ids=ids(stim.ResponseData.Values(2,:)==1);
    end
    line([round(nanmean(stim.ResponseData.Values(1,ids))) round(nanmean(stim.ResponseData.Values(1,ids)))],[0 20],'Color',[0 0 0.8]*i./4,'LineWidth',2);
end
xlabel('Reaction time [ms]');
ylabel('Number of trials [n]');
legend data MeanOfAll Control Famous Familiar Self 0.22 0.3 0.45 0.55

%% ERPs
for ch=0:8:56
    % for ch=0:1
    figure;
    c=0;
    for chf=1+ch:8+ch
        c=c+1;
        cc=0;
        %         indx=1:size(ALLEEG(ds).data,3); %all
        if cats==1
            for counte=[1 2 4 3]
                colors=[0 0.45 0.7 1];
                cc=cc+1;
                indx{cc}=[stim.stimTrain.imageCategory]==counte;
                if only_answered_trials==1
                    tmp=indx{cc};
                    tmp(isnan(stim.ResponseData.Values(2,:)))=0;
                    indx{cc}=tmp;
                    if only_correct_trials==1
                        tmp=indx{cc};
                        tmp(stim.ResponseData.Values(2,:)==0)=0;
                        indx{cc}=tmp;
                    end
                end
                subplot(2,4,c);
                plot([epoch_span(1):epoch_span(2)],squeeze(mean(ALLEEG(ds).data(chf,:,indx{cc}),3)),'Color',[0.8 0.2 1]*colors(cc),'LineWidth',1.5);
                hold on;
            end
        else
            for counte=[0.22 0.3 0.45 0.55]
                colors=[0 0.45 0.7 1];
                cc=cc+1;
                indx{cc}=[stim.stimTrain.imageNoise]==counte;
                if only_answered_trials==1
                    tmp=indx{cc};
                    tmp(isnan(stim.ResponseData.Values(2,:)))=0;
                    indx{cc}=tmp;
                    if only_correct_trials==1
                        tmp=indx{cc};
                        tmp(stim.ResponseData.Values(2,:)==0)=0;
                        indx{cc}=tmp;
                    end
                end
                subplot(2,4,c);
                plot([epoch_span(1):epoch_span(2)],squeeze(mean(ALLEEG(ds).data(chf,:,indx{cc}),3)),'Color',[1 0.2 0.8]*colors(cc),'LineWidth',1.5);
                hold on;
            end
        end
        line([epoch_span(1) epoch_span(2)],[0 0]);
        line([0 0],[-30 30]);
        xlim([epoch_span(1) epoch_span(2)]);
%         ylim([-15 15]);
        xlabel('Time [ms]');
        ylabel('Amplitude [uv]');
        title(EEG.chanlocs(chf).labels)
        hold off;
    end
    if cats==1
        legend ('Control','Famous','Familiar','Self','Location','southeast');
    else
        legend ('Noise = 0.22','Noise = 0.30','Noise = 0.45','Noise = 0.55','Location','southeast');
    end
end

%% Whole brain ERPs
chf=[1:64];
cc=0;
for counte1=[0.22 0.3 0.45 0.55]
    cc=cc+1;
    dd=0;
    for counte2=[1 2 4 3]
        dd=dd+1;
        indx{cc,dd}=(([stim.stimTrain.imageNoise]==counte1) & ([stim.stimTrain.imageCategory]==counte2));
        if only_answered_trials==1
            tmp=indx{cc,dd};
            tmp(isnan(stim.ResponseData.Values(2,:)))=0;
            indx{cc,dd}=tmp;
            if only_correct_trials==1
                tmp=indx{cc,dd};
                tmp(stim.ResponseData.Values(2,:)==0)=0;
                indx{cc,dd}=tmp;
            end
        end
        subplot(2,2,cc);
        plot([epoch_span(1):epoch_span(2)],squeeze(mean(mean(ALLEEG(ds).data(chf,:,indx{cc,dd}),3),1)),'LineWidth',1.5);        
        title(num2str(counte1));
        hold on;
    end
    line([epoch_span(1) epoch_span(2)],[0 0]);
    line([0 0],[-30 30]);
    xlim([epoch_span(1) epoch_span(2)]);
%     ylim([-15 15]);
    xlabel('Time [ms]');
    ylabel('Amplitude [uv]');
    hold off;
end
if cats==1
    legend ('Control','Famous','Familiar','Self','Location','southeast');
else
    legend ('Noise = 0.22','Noise = 0.30','Noise = 0.45','Noise = 0.55','Location','southeast');
end




