clc;
close all;
clear all;

region=1; %occ=1; parti=2; wholehead=3;
stim_resp=1; %stim_aligned=1, resp_aligned

for trials=[3 2 1]  %all=3; cor=2; incor=1;
    clearvars Means Vars Cors
    if trials==3
        if stim_resp==1
            load(['st_ParametersV2_All_trials_region_',num2str(region),'.mat']);
        else
            load(['rp_ParametersV2_All_trials_region_',num2str(region),'.mat']);
        end
    elseif trials==2
        if stim_resp==1
            load(['st_ParametersV2_Cor_trials_region_',num2str(region),'.mat']);
        else
            load(['rp_ParametersV2_Cor_trials_region_',num2str(region),'.mat']);
        end
    elseif trials==1
        if stim_resp==1
            load(['st_ParametersV2_Incor_trials_region_',num2str(region),'.mat']);
        else
            load(['rp_ParametersV2_Incor_trials_region_',num2str(region),'.mat']);
        end
    end
    
    smoothing=3;
    t=1:195;
    steps=10;
    
    if stim_resp==1
        for time=1:195
            Means(:,:,time,:)=Means(:,:,time,:)-repmat(nanmean(nanmean(Means(:,:,1:45,:),3),2),[1 4 1 1]);
        end
    end
    
    FanoFactor=Vars./abs(Means);
    for ch=1:size(FanoFactor,1)
        for cond=1:size(FanoFactor,2)
            for subj=[1:16 20 21]
                FanoFactor(ch,cond,isoutlier(FanoFactor(ch,cond,:,subj)),subj)=nan;
            end
        end
    end
    
    
    
    if stim_resp==1
        span=[-500:steps:1440]+25;
    else
        span=[-1450:steps:490]-25;
    end
    
    separatedByCats=0;
    if separatedByCats==1
        figure;
        for cat=1:4
                            plot(span,squeeze(smooth(nanmean(nanmean(FanoFactor(:,cat,t,:),1),4),smoothing)),'linewidth',3)
%             plot(span,squeeze(smooth(nanmean(nanmean(nanmean(Cors(:,:,cat,:,:),5),1),2),smoothing)),'linewidth',3)
            hold on;
        end
    else
                    plot(span,squeeze(smooth(nanmean(nanmean(nanmean(FanoFactor(:,:,t,:),2),1),4),smoothing)),'linewidth',3)
%         plot(span,squeeze(smooth(nanmean(nanmean(nanmean(nanmean(Cors,5),1),2),3),smoothing)),'linewidth',3)
        hold on;
    end
end




line([0 0],[0.7 0.9])
grid on;

if separatedByCats==1
    if stim_resp==1
        legend ('Control','Famous','Pers Familiar','Self','Location','northeast');
    else
        legend ('Control','Famous','Pers Familiar','Self','Location','northwest');
    end
else
    if stim_resp==1
        legend ('All trials','Cor trials','Incor trials','Location','northeast');
    else
        legend ('All trials','Cor trials','Incor trials','Location','northwest');
    end

end
xlabel('Time [ms]')
ylabel('Fano factor [Var/Mean]')
% ylabel('Corr [r]')
