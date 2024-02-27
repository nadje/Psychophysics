function SoEThresh_Detection
    % Behavioral SoEThresh experiment (experiment 1c of TEMPLATES project)
    % Financed by TEMPLATES project (PSI2014- 52573-P) and RYC-2013-12577 (MINECO)
    % PI: Iria SanMiguel

    % Experiment run by Raquel Aparicio (Bachelor thesis) and Nadia Paraskevoudi
    % Start July 2018 - End January 2019

    % This script aggregates all the 40 blocks from the detection task for each
    % individual subject and fits the psychometric functions for active and
    % passive conditions using the Palamedes toolbox. A psychometric curve is
    % fitted to the obtained datapoints using a Cumulative Normal function.
    % Fitting is done employing the PALAMEDES toolbox
    % (http://www.palamedestoolbox.org/), using the "PAL_PFML_Fit.m" function.
    % See adjustable parameters in the script.

    % ANALYSIS STEPS:

    % The script loads all the blocks from the results folder and merges them
    % in one mat structure (results). It rearranges the results to create two
    % separate matrices, one for active and one for passive, and for each
    % matrix: It gets rid of trials where no response was given. It gets a
    % vector with the unique intensity-values It creates a new matrix for each
    % condition to store the data that will be used for fitting the function.
    % This matrix consists of 4 columns with the following information:
    % Column1= intensity values, Column2 = Sum of the �Sound in Interval 2�
    % responses for each intensity, Column3 = total number of times that this
    % intensity was presented, Column4= absolute number of correct responses
    % for this intensity, Column5= Percent correct.

    % "results_no_miss" matrix (stored as e.g. "99_Detection_Passive_Raw.mat"
    % and "99_Detection_Active_Raw.mat", respectively) To create this matrix,
    % we take each results matrix saved for each block and merge the 25
    % matrices for each participant into one matrix for the active and passive
    % conditions separately. Thus, the new matrix is a merged version of all
    % the blocks. The new merged matrix does not contain trials without
    % response (i.e., response == 3), since these trials will not be used for
    % the analysis

    % "results_forfit" matrix is created separately for the active and passive
    % trials(stored in e.g."99_Detection_Passive.mat" and
    % "99_Detection_Active.mat", respectively) Contains results in format ready
    % for psychometric function fit later in R. Each row is one
    % sound intensity presented at least once, and in columns: 1) Intensity in
    % dB 2) Total number of times this intensity has been presented 3) Number
    % of correct responses for each intensity level 4) Percent correct for each
    % intensity level

    % OUTPUT: Results matrix (raw results) and results_forfit matrix
    % (results processed for psychometric fit function) - all results are
    % automatically stored as files in the /res folder, named with the subject
    % number

    % AFTER THE PER-SUBJECT ANALYSIS
    % Load data matrices for each subject and merge them into a big file with
    % participants ID. This data will be later introduced in R for the fitting
    % and simulation procedures.

    % Script/Comments Updated by Nadia - December 3rd 2019

    analysis_scripts = 'D:\SoEThresh_scripts';
    results_path = 'D:\SoEThresh_scripts\Raw_Data\Detection_Task';
    analysis_path = 'D:\SoEThresh_scripts\Analysis';
    all_res = 'D:\SoEThresh_scripts\Analysis';

    nBlocks = 25; % Define the number of blocks for each participant

    % Load and aggregate result matrices from the two blocks
    for subject = 1:length(subjects)
        fprintf(anal_logfile, '\r\n%s', ['Analyzing subject ' num2str(subjects(subject))]); % update analysis log
        results_all = [];

        for iBlock = 1:nBlocks
            load([results_path '/' num2str(subjects(subject), '%.02d') '_' num2str(iBlock, '%.2d')]);
            results_all = [results_all; matlog];
        end

        results = results_all;

        %% CREATE TWO SEPARATE MATRICES FOR ACTIVE AND PASSIVE TRIALS
        %----------------------------------------------------------------------
        %               Rearrange results format for Palamedes
        %----------------------------------------------------------------------

        % results so far are stored in the matrix "matlog"
        % we need to sort and calculate the percent correct for each
        % intensity and for each task type (active and passive)
        results(results(:, 23) == 0, :) = []
        %%%%%%% ACTIVE TRIALS - RESULTS %%%%%%%%
        results_nomiss_active = results;
        results_nomiss_active(results_nomiss_active(:, 4) == 3, :) = []; % get rid of iterations without response (i.e., response == 3)
        results_nomiss_active(results_nomiss_active(:, 1) == 0, :) = []; % get rid of passive trials for this file - we only need the active trials
        int_values = unique(results_nomiss_active(:, 3)); % get vector of intensities delivered at least once, sorted
        results_forfit_active = zeros((length(int_values)), 4); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns
        results_forfit_active(:, 1) = int_values; % first column is intensity values

        for iInt = 1:length(int_values) %for each one of the unique intensity values presented
            results_forfit_active(iInt, 2) = sum(results_nomiss_active(:, 3) == int_values(iInt)); % second column is total number of times this intensity has been presented
            results_forfit_active(iInt, 3) = sum(results_nomiss_active(results_nomiss_active(:, 3) == int_values(iInt), 5)); %third column is number of correct responses for each intensity value
            results_forfit_active(iInt, 4) = (100 * results_forfit_active(iInt, 3)) / results_forfit_active(iInt, 2); %Percent correct
        end

        %%%%%%% PASSIVE TRIALS - RESULTS %%%%%%%%
        results_nomiss_passive = results;
        results_nomiss_passive(results_nomiss_passive(:, 4) == 3, :) = []; % get rid of iterations without response (i.e., response ==3)
        results_nomiss_passive(results_nomiss_passive(:, 1) == 1, :) = []; % get rid of active trials for this file - we only need the passive trials

        int_values = unique(results_nomiss_passive(:, 3)); % get vector of intensities delivered at least once, sorted
        results_forfit_passive = zeros(length(int_values), 3); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns

        results_forfit_passive(:, 1) = int_values; % first column is intensity values

        for ip = 1:length(int_values)
            results_forfit_passive(ip, 2) = sum(results_nomiss_passive(:, 3) == int_values(ip)); % second column is total number of times this intensity has been presented
            results_forfit_passive(ip, 3) = sum(results_nomiss_passive(results_nomiss_passive(:, 3) == int_values(ip), 5)); %third column is number of correct responses for each intensity value
            results_forfit_passive(ip, 4) = (100 * results_forfit_passive(ip, 3)) / results_forfit_passive(ip, 2); %Percent correct
        end

        %---------------------------------------------------------------------
        %              SAVE
        %----------------------------------------------------------------------
        cd(analysis_path)

        % save results_forfit_active matrix
        eval(['save ' num2str(subjects(subject), '%02d') '_Detection_Active results_forfit_active']);
        % save results_nomiss_passive (with timing info)
        eval(['save ' num2str(subjects(subject), '%02d') '_Detection_Active_Raw results_nomiss_active']);

        % save results_forfit_passive matrix
        eval(['save ' num2str(subjects(subject), '%02d') '_Detection_Passive results_forfit_passive']);
        % save results_nomiss_passive (with timing info)
        eval(['save ' num2str(subjects(subject), '%02d') '_Detection_Passive_Raw results_nomiss_passive']);
        cd(analysis_scripts)

    end

    %% WHOLE SAMPLE ANALYSIS

    % In the array "subjects" only the final sample needs to be included
    % Participants 2, 19, 25 excluded. Final Sample N = 28
    subjects = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 31];
    cd(analysis_scripts);

    res_active = [];
    res_passive = [];
    res_activeRaw = [];
    res_passiveRaw = [];
    % (active and passive)
    for iSubject = 1:length(subjects)
        % First load individuals' thresholds
        load([analysis_path '/' num2str(subjects(iSubject), '%.2d') '_Detection_Active'])
        load([analysis_path '/' num2str(subjects(iSubject), '%.2d') '_Detection_Passive'])
        load([analysis_path '/' num2str(subjects(iSubject), '%.2d') '_Detection_Active_Raw'])
        load([analysis_path '/' num2str(subjects(iSubject), '%.2d') '_Detection_Passive_Raw'])
        % Then to new collapsed matrices
        results_nomiss_active = results_nomiss_active(:, [1:5]); %keep only columns 1 to 5
        results_nomiss_passive = results_nomiss_passive(:, [1:5]);
        results_nomiss_active(:, 6) = repmat(subjects(iSubject), 1, length(results_nomiss_active)) % add participants' ID
        results_nomiss_passive(:, 6) = repmat(subjects(iSubject), 1, length(results_nomiss_passive))
        res_activeRaw = [res_activeRaw; results_nomiss_active];
        res_passiveRaw = [res_passiveRaw; results_nomiss_passive];

        results_forfit_active(:, 5) = repmat(subjects(iSubject), 1, 8); % add Participants' ID
        results_forfit_passive(:, 5) = repmat(subjects(iSubject), 1, 8); % add Participants' ID
        res_active = [res_active; results_forfit_active];
        res_passive = [res_passive; results_forfit_passive];
    end

    cd(analysis_path)
    %    Save
    eval(['save ' 'res_active res_active']);
    eval(['save ' 'res_passive res_passive']);
    eval(['save ' 'res_activeRaw res_activeRaw']);
    eval(['save ' 'res_passiveRaw res_passiveRaw']);

    fclose(anal_logfile);
end
