function SoEThresh_Detection_DPrimeAnalysis
% based on that https://www.birmingham.ac.uk/Documents/college-les/psych/vision-laboratory/sdtintro.pdf
% with all intensities included

% The d' prime is calculated as follows:

%--------------------------- Trial order ----------------------------------
%              Signal-Noise (1st int) |   Noise - Signal (2nd int)
%-------------------------------------|------------------------------------
% R | 1st            Hit              |             False alarm
% E |---------------------------------|------------------------------------
% S |2nd             1- Hit           |             1 - FA
% P |----------------------------------------------------------------------

% Saves Hit Rate, False Alarm Rate, Percent Correct, d'prime, and criterion
% per participant

% In the final "Whole sample analysis", it creates a big matrix including
% all participants for d', criterion, percent correct, separately and saves
% it in csv. format to be later introduced in R. 

% Updated by Nadia - December 3rd 2019
% -----------------------------------------------------

% DEFINE PATHS
analysis_scripts = 'D:\SoEThresh_scripts';
results_path = 'D:\SoEThresh_scripts\Raw_Data\Detection_Task';
analysis_path = 'D:\SoEThresh_scripts\Analysis';
all_res = 'D:\SoEThresh_scripts\Analysis';
anal_logfile = [analysis_path '\analysis_log.txt'];


cd(analysis_scripts);

% Include final sample
% Participants 2,19, 25 exluded
subjects_new = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 26,27,28,29,30,31];

for subject = 1:length(subjects_new)
    load([analysis_path '/'  num2str(subjects_new(subject),'%02d') '_Detection_Active_Raw']);
    load([analysis_path '/' num2str(subjects_new(subject),'%02d') '_Detection_Passive_Raw']);
    load([analysis_path '/' num2str(subjects_new(subject),'%02d') '_Detection_Passive']);
    load([analysis_path '/' num2str(subjects_new(subject),'%02d') '_Detection_Active']);
    
    TrialOrder = []; % 1 = Signal in int 1 and 2 = signal in Int 2
    iResp =[]; % 1 = Response1, 2=Response 2
    
    hit = false;
    FA=false;
    c = 1; %counter
    
    % Hit: when sound was in 1 and pressed 1 (cond =0 and resp=0)  (button 1)
    % Miss  when sound was in 1 and pressed 2 (cond =0 and resp=1)
    
    % False Alarmns    when sound was in 2 and pressed 1 (cond =1 and resp=0)
    % Correct Rejections  when sound was in 2 and pressed 2 (cond =0 and resp=0)
    n_active = length(results_nomiss_active(:,1));
    n_passive = length(results_nomiss_passive(:,1));
    dprime_active = zeros(n_active,5);
    dprime_passive = zeros(n_passive,5);
    
    %%%%% ACTIVE TRIALS %%%%%%%%%%%%%%%
    for c = 1:n_active
        db = results_nomiss_active(c,3);
        cond = results_nomiss_active(c,2);
        resp = results_nomiss_active(c,4);
        TrialOrder = 0;
        iResp = 0;
        hit = false;
        
        FA=false;
        
        if cond == 0 && resp == 0 %hits
            TrialOrder = 0;
            iResp = 0;
            hit = true;
        elseif cond == 1 && resp == 0 % false alarms
            TrialOrder = 1;
            iResp = 0;
            FA =true;
        end
        
        dprime_active(c, 1) = db;
        dprime_active(c, 2) = TrialOrder;
        dprime_active(c, 3) = iResp;
        dprime_active(c, 4) = hit
        dprime_active(c, 5) = FA;
    end
     nHits = sum(dprime_active(:,4) ==1)
    nMisses = length(dprime_active(:,4)) - nHits
    nFA = sum(dprime_active(:,5) ==1);
  nCR = length(dprime_active(:,5)) - nFA
  
  active_for_R = [nHits, nMisses, nFA, nCR ,length(dprime_active(:,5))] % save for R (June 2020)
    % Calculate Hit Rate, False alarm Rate, Percent Correct, d', and
    % criterion for Active trials
    HR_A =((100*sum(dprime_active(:,4)))/sum(results_forfit_active(:,2)))/100
    FAR_A = ((100*sum(dprime_active(:,5)))/sum(results_forfit_active(:,2)))/100
    Pc_A = (HR_A+(1-FAR_A))/2
    dpA = norminv(HR_A)-norminv(FAR_A) %d'prime 
    c_A = -0.5*(norminv(HR_A)+ norminv(FAR_A)); %criterion
    
    %%%%% PASSIVE TRIALS %%%%%%%%%%%%%%%
    for c = 1:n_passive
        db = results_nomiss_passive(c,3);
        cond = results_nomiss_passive(c,2);
        resp = results_nomiss_passive(c,4);
        TrialOrder = 0;
        iResp = 0;
        hit = false;
        
        FA=false;
        
        if cond == 0 && resp == 0 %hits
            TrialOrder = 0;
            iResp = 0;
            hit = true;
        elseif cond == 1 && resp == 0 % false alarms
            TrialOrder = 1;
            iResp = 0;
            FA =true;
        end
        dprime_passive(c, 1) = db;
        dprime_passive(c, 2) = TrialOrder;
        dprime_passive(c, 3) = iResp;
        dprime_passive(c, 4) = hit
        dprime_passive(c, 5) = FA;
    end
    
      nHits = sum(dprime_passive(:,4) ==1)
    nMisses = length(dprime_passive(:,4)) - nHits
    nFA = sum(dprime_passive(:,5) ==1);
  nCR = length(dprime_passive(:,5)) - nFA
  
  passive_for_R = [nHits, nMisses, nFA, nCR, length(dprime_passive(:,5))] % save for R (June 2020)
    % Calculate Hit Rate, False alarm Rate, Percent Correct, d', and
    % criterion for Passive trials
    HR_P = ((100*sum(dprime_passive(:,4)))/sum(results_forfit_passive(:,2)))/100
    FAR_P = ((100*sum(dprime_passive(:,5)))/sum(results_forfit_passive(:,2)))/100
    Pc_P = (HR_P+(1-FAR_P))/2
    dpP = norminv(HR_P)-norminv(FAR_P)
    c_P = -0.5*(norminv(HR_P)+ norminv(FAR_P))
    
    cd(analysis_path);
    
    %Save Hit Rate & False alarm rate
    eval(['save ' num2str(subjects_new(subject),'%02d') '_active_for_R active_for_R']);
        eval(['save ' num2str(subjects_new(subject),'%02d') '_passive_for_R passive_for_R']);

    
    eval(['save ' num2str(subjects_new(subject),'%02d') '_HR HR_A HR_P']);
    eval(['save ' num2str(subjects_new(subject),'%02d') '_FA FAR_A FAR_P']);
    %Save d prime, criterion, and percent correct
    eval(['save ' num2str(subjects_new(subject),'%02d') '_dprimes dpA dpP']);
    eval(['save ' num2str(subjects_new(subject),'%02d') '_criterion c_A c_P'])
    eval(['save ' num2str(subjects_new(subject),'%02d') '_pCorrect Pc_A Pc_P'])
   
    end
    
    

% WHOLE SAMPLE ANALYSIS

% Create a big file for each measure to be introduced later in R. 
dps =[];
crits =[];
pcs =[];
activeSDT = []; passiveSDT=[];
for iSubject = 1:length(subjects_new)
    
    load([analysis_path '/' num2str(subjects_new(iSubject),'%.2d') '_dprimes']); % results-for-fit passive
    load([analysis_path '/' num2str(subjects_new(iSubject),'%.2d') '_criterion']); % results-for-fit active
    load([analysis_path '/' num2str(subjects_new(iSubject),'%.2d') '_pCorrect']); % results-for-fit passive
        load([analysis_path '/' num2str(subjects_new(iSubject),'%.2d') '_active_for_R']); % results-for-fit passive
  load([analysis_path '/' num2str(subjects_new(iSubject),'%.2d') '_passive_for_R']); % results-for-fit passive
active_for_R(end+1) = subjects_new(iSubject);
activeSDT = [activeSDT; active_for_R];

passive_for_R(end+1) = subjects_new(iSubject);
passiveSDT = [passiveSDT; passive_for_R];


    dps = [dps; dpA dpP]
    crits =[crits; c_A c_P]
    pcs =[pcs; Pc_A Pc_P];
    
end

% SAVE the merged files in csv format to be used later in R
cd(all_res);
save passiveSDT
save activeSDT
csvwrite('dprimes_for_R.csv', dps);
csvwrite('criterion_for_R.csv', crits);
csvwrite('pCorrect_for_R.csv', pcs);
end


