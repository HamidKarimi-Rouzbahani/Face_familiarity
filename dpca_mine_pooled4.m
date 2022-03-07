clc;
clear all;
% close all;

%% Settings
% cats=0; % categories=1 and coh=0
only_answered_trials=1; % all trials=0

% stim_resp=1; % for stim=1, for resp=2
baseline_correction=1;

%% ERPs
subject={'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
sessions=[1 2];
ifSimultaneousRecording = false;

for region=[1]
    for stim_resp=[1]
        clearvars -except region stim_resp ifSimultaneousRecording sessions subject baseline_correction only_answered_trials
        firingRatesPooled=nan*ones(1,4,2,2000,240);
        firingRatesAveragePooled=nan*ones(1,4,2,2000);
        trialNumPooled=zeros(1,4,2);
        
        for subj=[1:16 20 21]
            if region==3
                firingRates=nan*ones(62,4,2,2000,240);
                trialNum=zeros(62,4,2);
            else
                firingRates=nan*ones(14,4,2,2000,240);
                trialNum=zeros(14,4,2);
            end
            for correct_incorrect_trials=[1:2]
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
                        indx{cc,session}=[Exp_data{1,session}.stim.stimTrain.imageCategory]==counte;
                        if only_answered_trials==1
                            tmp=indx{cc,session};
                            tmp(isnan(Exp_data{1,session}.stim.ResponseData.Values(2,:)))=0;
                            indx{cc,session}=tmp;
                            if correct_incorrect_trials==1
                                tmp=indx{cc,session};
                                tmp(Exp_data{1,session}.stim.ResponseData.Values(2,:)==0)=0;
                                indx{cc,session}=tmp;
                            elseif correct_incorrect_trials==2
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
                
                
                if baseline_correction && stim_resp==1
                    for i=1:4
                        signals{1,i}=signals{1,i}-repmat(nanmean(nanmean(signals{1,i}(:,baseline_span,:),3),2),[1 2000 size(signals{1,i},3)]);
                    end
                end
                % if cats==1
                %     legend ('Control','Famous','Familiar','Self','Location','southeast');
                % else
                %     legend ('Coherence = 0.22','Coherence = 0.30','Coherence = 0.45','Coherence = 0.55','Location','southeast');
                % end
                
                if region==1
                    for i=1:4
                        signals{1,i}([1:12 13 14:23 25:27 32 33:49 51:52 58:61],:,:)=[];
                    end
                elseif region==2
                    for i=1:4
                        signals{1,i}([1:5 7:9 12 13 14:15 17:20 23:25 27:31 32 33:37 40:41 43:44 47:50 53:64],:,:)=[];
                    end
                elseif region==3
                    for i=1:4
                        signals{1,i}([13 32],:,:)=[];
                    end
                end
                
%                 % Adding all familiar categories to cond #2
%                 signals{1,2}(:,:,end+1:end+size(signals{1,3},3))=signals{1,3};
%                 signals{1,2}(:,:,end+1:end+size(signals{1,4},3))=signals{1,4};
%                 signals{1,3}=[];
%                 signals{1,4}=[];
                size_signals=4;
%                 for j=1:length(signals)
%                     if ~isempty(signals{1,j})
%                         size_signals=size_signals+1;
%                     end
%                 end
                for j=1:size_signals
                    firingRates(:,j,correct_incorrect_trials,1:2000,1:size(signals{1,j},3))=signals{1,j};
                    trialNum(:,j,correct_incorrect_trials)=size(signals{1,j},3);
                end
                ccc
                clearvars -except firingRatesPooled trialNumPooled firingRatesAveragePooled region ifSimultaneousRecording trialNum firingRates test_on_errors sizes_signals accuracy Exp_data baseline_correction EEG_signals subject sessions steps baseline_span signals2 signals subj cats stim_resp only_answered_trials only_correct_trials only_incorrect_trials only_answered_trials2 only_correct_trials2 only_incorrect_trials2
            end
            firingRatesAverage = nanmean(firingRates, 5);
            
            firingRatesPooled=cat(1, firingRatesPooled,firingRates);
            trialNumPooled=cat(1, trialNumPooled,trialNum);
            firingRatesAveragePooled=cat(1,firingRatesAveragePooled,firingRatesAverage);
        end
        
        firingRates=firingRatesPooled(2:end,:,:,:,:);
        firingRatesAverage=firingRatesAveragePooled(2:end,:,:,:);
        trialNum=trialNumPooled(2:end,:,:);
        
        N = size(firingRates,1);    % number of neurons
        S = size(firingRates,2);    % number of stimuli
        D = size(firingRates,3);    % number of decisions
        T = size(firingRates,4);    % number of time points
        E = size(firingRates,5);    % maximal number of trial repetitions
        %% Define parameter grouping
        
        % *** Don't change this if you don't know what you are doing! ***
        % firingRates array has [N S D T E] size; herewe ignore the 1st dimension
        % (neurons), i.e. we have the following parameters:
        %    1 - stimulus
        %    2 - decision
        %    3 - time
        % There are three pairwise interactions:
        %    [1 3] - stimulus/time interaction
        %    [2 3] - decision/time interaction
        %    [1 2] - stimulus/decision interaction
        % And one three-way interaction:
        %    [1 2 3] - rest
        % As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:
        
        combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
        margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
        
        % For two parameters (e.g. stimulus and time, but no decision), we would have
        % firingRates array of [N S T E] size (one dimension less, and only the following
        % possible marginalizations:
        %    1 - stimulus
        %    2 - time
        %    [1 2] - stimulus/time interaction
        % They could be grouped as follows:
        %    combinedParams = {{1, [1 2]}, {2}};
        
        % Time events of interest (e.g. stimulus onset/offset, cues etc.)
        % They are marked on the plots with vertical lines
        if stim_resp==1
            time = ([-0.5:0.001:1.499]);
            timeEvents = time(round(length(time)/4));
        else
            time = ([-1.499:0.001:0.5]);
            timeEvents = time(round(length(time)*3/4));
        end
        % check consistency between trialNum and firingRates
        for n = 1:size(firingRates,1)
            for s = 1:size(firingRates,2)
                for d = 1:size(firingRates,3)
                    assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
                end
            end
        end
        %% Step 3: dPCA without regularization and ignoring noise covariance
        
        % This is the core function.
        % W is the decoder, V is the encoder (ordered by explained variance),
        % whichMarg is an array that tells you which component comes from which
        % marginalization
        
        tic
        NumberOfEigens=min(20,N);
        
        [W,V,whichMarg] = dpca(firingRatesAverage, NumberOfEigens, ...
            'combinedParams', combinedParams);
        toc
        [region stim_resp]
        save(['DPCA_Cats_region_',num2str(region),'_stim_resp_',num2str(stim_resp),'.mat'],'W','V','whichMarg');
        
%         %% plotting
%         load(['DPCA_Cats_region_',num2str(region),'_stim_resp_',num2str(stim_resp),'.mat']);
%         
%         explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
%             'combinedParams', combinedParams);
%         
%         dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%             'explainedVar', explVar, ...
%             'marginalizationNames', margNames, ...
%             'marginalizationColours', margColours, ...
%             'whichMarg', whichMarg,                 ...
%             'time', time,                        ...
%             'timeEvents', timeEvents,               ...
%             'timeMarginalization', 3, ...
%             'legendSubplot', 16);
    end
end
