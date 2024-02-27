function SoEThresh_Discrimination
    % Behavioral SoEThresh experiment (experiment 1c of TEMPLATES project)
    % Financed by TEMPLATES project (PSI2014- 52573-P) and RYC-2013-12577 (MINECO)
    % PI: Iria SanMiguel

    % Experiment run by Raquel Aparicio (Bachelor thesis) and Nadia Paraskevoudi
    % Start July 2018 - End January 2019

    % This script aggregates all the 25 blocks from the discrimination task for
    % each individual subject. The fitting and simulation/bootstrapping is done
    % later in R.

    % Analysis is based on the following steps:

    % A) The script loads all the blocks from the results folder and merges
    % them in one mat structure (results). It rearranges the results to create
    % 4 separate matrices, one for each condition. Then, for each matrix: It
    % gets rid of trials where no response was given, as well sas of trials
    % where participants pressed the control-X button to indicate that they did
    % not listen to the sound.

    % B)It gets a vector with the unique intensity-values for each condition It
    % creates a new matrix for each condition to store the data that will be
    % used for fitting the function. Instead of using the absolute values, we
    % calculated the relative values [-3:3] to be able to compare across supra-
    % and near-threhold conditions. This matrix consists of 4 columns with the
    % following information: Column1= absolute intensity values, Column2 = Sum
    % of the �Second sound louder� responses for each comparison intensity,
    % Column3 = total number of times that this intensity was presented,
    % Column4= absolute number of correct responses for this intensity,
    % Column5= Percent of �second sound louder� responses, Column5 = relative
    % intensity values [-3:3]

    % C) Calculate Percent correct for each subject/condition for when the two
    % intensities where the same. The following part of the analysis is
    % conducted so as to directly compare our results with the findings
    % reported by Reznik et al., 2015 Reznik et al. (2015) employed the exact
    % same paradigm, with the only difference that unbeknownst to the
    % participants, the two sounds (standard and comparison) was presented at
    % the same intensity. To directly compare with their study, we, thus,
    % rearrange the results by taking into account only the trials where the
    % comparison tone was presented at the same intensity as the standard one
    % (i.e., relative difference/ relative intensity value = 0). The rest of
    % the analysis is similar to the analysis described above.

    % OUTPUT: Saves 1) raw data for each condition separately, 2) matrix for
    % fitting, 3) matrix with results rearranged so as to compare with Reznik
    % et al. (2015) results.

    % AFTER THE PER-SUBJECT ANALYSIS
    % Load data matrices for each subject and merge them into a big file with
    % participants ID. This data will be later introduced in R for the fitting
    % and simulation procedures.

    % Script/Comments Updated by Nadia - December 3rd 2019

    % ------------------------------------------------------------------------

    % Define Paths
    experiment_path = 'D:\SoEThresh_scripts';
    analysis_scripts = 'D:\SoEThresh_scripts';
    results_path = 'D:\SoEThresh_scripts\Raw_Data\Discrimination_Task';
    analysis_path = 'D:\SoEThresh_scripts\Analysis';
    all_res = 'D:\SoEThresh_scripts\Analysis';

    anal_logfile = [analysis_path '/analysis_log.txt'];
    anal_logfile = fopen(anal_logfile, 'a');
    fprintf(anal_logfile, '\r\n%s', '******************************************');
    fprintf(anal_logfile, '\r\n%s', datestr(now)); % write date and time
    fprintf(anal_logfile, '\r\n%s', 'SoEThresh_Discrimination'); % Indicate which script
    fprintf(anal_logfile, '\r\n%s', ['Mat files will be stored in ' num2str(analysis_path)]); %
    nBlocks = 25; % number of blocks to be merged per participant
    subjects = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 31];

    for subject = 1:length(subjects)

        fprintf(anal_logfile, '\r\n%s', '--------------------------'); %
        fprintf(anal_logfile, '\r\n%s', ['Analyzing subject ' num2str(subjects(subject))]); % update analysis log
        fprintf(anal_logfile, '\r\n%s', '--------------------------'); %
        % Load and aggregate result matrices from the two blocks
        results_all = [];

        for iBlock = 1:nBlocks
            load([results_path '/' num2str(subjects(subject), '%.2d') '_' num2str(iBlock, '%.2d')]);
            results_all = [results_all; discr_log];
        end

        results = results_all;
        %% CREATE 4 SEPARATE MATRICES FOR EACH CONDITION
        %
        % Rearrange results format for Palamedes:
        % Results so far are stored in the matrix "matlog"
        % We need to sort and calculate the % of 2nd Sound louder responses for
        % each comparison intensity and for each condition.

        % -----------------------------------------------------------------------
        fprintf(anal_logfile, '\r\n%s', ' ');
        fprintf(anal_logfile, '\r\n%s', 'Misses and X-button presses');
        fprintf(anal_logfile, '\r\n%s', ' ');
        fprintf(anal_logfile, '\r\n%s', ['Misses ' num2str(sum(results(:, 4) == 3))]);
        fprintf(anal_logfile, '\r\n%s', ['X-button ' num2str(sum(results(:, 4) == 5))]);
        fprintf(anal_logfile, '\r\n%s', ' ');

        % ACTIVE SUPRATHRESHOLD TRIALS - RESULTS

        results_nomiss_as = results;

        % get rid of the other conditions (1,2,4)
        results_nomiss_as(results_nomiss_as(:, 1) == 1, :) = [];
        results_nomiss_as(results_nomiss_as(:, 1) == 2, :) = [];
        results_nomiss_as(results_nomiss_as(:, 1) == 4, :) = [];
        results_nomiss_as(results_nomiss_as(:, 1) == 0, :) = []; % get rid of a stupid zero-value array that insists on appearing
        % How many misses and X-responses for AS?
        tbd_as = sum(results_nomiss_as(:, 4) == 3) + sum(results_nomiss_as(:, 4) == 5) % all trials excluded
        tbd_x_as = sum(results_nomiss_as(:, 4) == 5); % only x button press for this condition
        tbd_miss_as = sum(results_nomiss_as(:, 4) == 3); % only missed responses for this condition

        % Now delete them from further analysis
        results_nomiss_as(results_nomiss_as(:, 4) == 3, :) = []; % get rid of iterations without response (i.e., response ==3)
        results_nomiss_as(results_nomiss_as(:, 4) == 5, :) = []; % get rid of the trials where the participant did not listen to the sound (i.e., response ==5)

        % Prepare matrix for the later fitting procedure
        int_values = unique(results_nomiss_as(:, 3)); % get vector of intensities delivered at least once, sorted

        results_forfit_as = zeros((length(int_values)), 6); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns
        results_forfit_as(:, 1) = int_values; % first column is intensity values

        for iInt = 1:length(int_values) %for each one of the unique intensity values presented
            results_forfit_as(iInt, 2) = sum(results_nomiss_as(results_nomiss_as(:, 3) == int_values(iInt), 4)); % second column is sum of the "2nd Sound is louder" responses for this intensity
            results_forfit_as(iInt, 3) = sum(results_nomiss_as(:, 3) == int_values(iInt)); % third column is total number of times this intensity has been presented
            results_forfit_as(iInt, 4) = (results_forfit_as(iInt, 2) ./ results_forfit_as(iInt, 3)) * 100; %percent of 2nd sound louder responses
        end

        %Create column with relative comparison values when compared with the standard intensity
        results_forfit_as(:, 5) = [-3:3];
        % -----------------------------------------------------------------------

        %  PASSIVE SUPRATHRESHOLD TRIALS - RESULTS

        results_nomiss_ps = results;

        % get rid of the rest of the conditions trials
        results_nomiss_ps(results_nomiss_ps(:, 1) == 2, :) = [];
        results_nomiss_ps(results_nomiss_ps(:, 1) == 3, :) = [];
        results_nomiss_ps(results_nomiss_ps(:, 1) == 4, :) = [];

        results_nomiss_ps(results_nomiss_ps(:, 3) == 0, :) = []; %get rid of this stupid empty line with zeros.
        % Find how many misses and discarded for PS

        tbd_ps = sum(results_nomiss_ps(:, 4) == 3) + sum(results_nomiss_ps(:, 4) == 5)
        tbd_x_ps = sum(results_nomiss_ps(:, 4) == 5); % only x button press for this condition
        tbd_miss_ps = sum(results_nomiss_ps(:, 4) == 3); % only missed responses for this condition
        % exclude these trials
        results_nomiss_ps(results_nomiss_ps(:, 4) == 3, :) = []; % get rid of iterations without response (i.e., response ==3)
        results_nomiss_ps(results_nomiss_ps(:, 4) == 5, :) = []; % get rid of the trials where the participant did not listen to the sound (i.e., response ==5)

        int_values = unique(results_nomiss_ps(:, 3)); % get vector of intensities delivered at least once, sorted

        results_forfit_ps = zeros(length(int_values), 4); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns

        results_forfit_ps(:, 1) = int_values; % first column is intensity values

        for ip = 1:length(int_values) %for each one of the unique intensity values presented

            results_forfit_ps(ip, 2) = sum(results_nomiss_ps(results_nomiss_ps(:, 3) == int_values(ip), 4)); % second column is sum of the "Sound in Interval 2" responses for this intensity
            results_forfit_ps(ip, 3) = sum(results_nomiss_ps(:, 3) == int_values(ip)); % third column is total number of times this intensity has been presented
            results_forfit_ps(ip, 4) = (results_forfit_ps(ip, 2) ./ results_forfit_ps(ip, 3)) * 100; %percent of 2nd sound louder responses
        end

        results_forfit_ps(:, 5) = [-3:3];
        % -----------------------------------------------------------------------
        %  ACTIVE NEARTHRESHOLD TRIALS - RESULTS

        results_nomiss_an = results;

        % get rid of the rest of the conditions trials
        results_nomiss_an(results_nomiss_an(:, 1) == 1, :) = [];
        results_nomiss_an(results_nomiss_an(:, 1) == 2, :) = [];
        results_nomiss_an(results_nomiss_an(:, 1) == 3, :) = [];

        results_nomiss_an(results_nomiss_an(:, 3) == 0, :) = []; %get rid of this stupid empty line with zeros.
        % Misses and X-button presses for AN
        tbd_an = sum(results_nomiss_an(:, 4) == 3) + sum(results_nomiss_an(:, 4) == 5)
        tbd_x_an = sum(results_nomiss_an(:, 4) == 5); % only x button press for this condition
        tbd_miss_an = sum(results_nomiss_an(:, 4) == 3); % only missed responses for this condition

        results_nomiss_an(results_nomiss_an(:, 4) == 3, :) = []; % get rid of iterations without response (i.e., response ==3)
        results_nomiss_an(results_nomiss_an(:, 4) == 5, :) = []; % get rid of the trials where the participant did not listen to the sound (i.e., response ==5)

        %Prepare matrix for fitting
        int_values = unique(results_nomiss_an(:, 3)); % get vector of intensities delivered at least once, sorted

        results_forfit_an = zeros(length(int_values), 4); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns

        results_forfit_an(:, 1) = int_values; % first column is intensity values

        for an = 1:length(int_values) %for each one of the unique intensity values presented

            results_forfit_an(an, 2) = sum(results_nomiss_an(results_nomiss_an(:, 3) == int_values(an), 4)); % second column is sum of the "Sound in Interval 2" responses for this intensity
            results_forfit_an(an, 3) = sum(results_nomiss_an(:, 3) == int_values(an)); % third column is total number of times this intensity has been presented
            results_forfit_an(an, 4) = (results_forfit_an(an, 2) ./ results_forfit_an(an, 3)) * 100; %Percent louder
        end

        %Create column with relative comparison values when compared with the standard intensity
        results_forfit_an(:, 5) = [-3:3];
        % -----------------------------------------------------------------------
        % PASSIVE NEARTHRESHOLD TRIALS - RESULTS

        results_nomiss_pn = results;
        % get rid of the rest of the conditions trials
        results_nomiss_pn(results_nomiss_pn(:, 1) == 1, :) = [];
        results_nomiss_pn(results_nomiss_pn(:, 1) == 3, :) = [];
        results_nomiss_pn(results_nomiss_pn(:, 1) == 4, :) = [];
        results_nomiss_pn(results_nomiss_pn(:, 1) == 0, :) = []; %get rid of this stupid empty line with zeros.
        % Misses and X-Button presses for PN
        tbd_pn = sum(results_nomiss_pn(:, 4) == 3) + sum(results_nomiss_pn(:, 4) == 5)
        tbd_x_pn = sum(results_nomiss_pn(:, 4) == 5); % only x button press for this condition
        tbd_miss_pn = sum(results_nomiss_pn(:, 4) == 3); % only missed responses for this condition

        results_nomiss_pn(results_nomiss_pn(:, 4) == 3, :) = []; % get rid of iterations without response (i.e., response ==3)
        results_nomiss_pn(results_nomiss_pn(:, 4) == 5, :) = [];

        % Prepare matrix for later fitting
        int_values = unique(results_nomiss_pn(:, 3)); % get vector of intensities delivered at least once, sorted
        results_forfit_pn = zeros(length(int_values), 4); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns
        results_forfit_pn(:, 1) = int_values; % first column is intensity values

        for pn = 1:length(int_values) %for each one of the unique intensity values presented

            results_forfit_pn(pn, 2) = sum(results_nomiss_pn(results_nomiss_pn(:, 3) == int_values(pn), 4)); % second column is sum of the "Sound 2 louder" responses for this intensity
            results_forfit_pn(pn, 3) = sum(results_nomiss_pn(:, 3) == int_values(pn)); % third column is total number of times this intensity has been presented
            results_forfit_pn(pn, 4) = (results_forfit_pn(pn, 2) ./ results_forfit_pn(pn, 3)) * 100; %percent of Sound 2 louder responses
        end

        %Create column with relative comparison values when compared with the standard intensity
        results_forfit_pn(:, 5) = [-3:3];

        % ----
        % Get a summary of the trials to be discarded per participant and
        % condition (Based on Reviewer 1 comment)
        number_of_trials_to_be_excluded(subject, :) = [tbd_as, tbd_ps, tbd_an, tbd_pn]
        number_of_trials_to_be_excluded_x(subject, :) = [tbd_x_as, tbd_x_ps, tbd_x_an, tbd_x_pn];
        number_of_trials_to_be_excluded_miss(subject, :) = [tbd_miss_as, tbd_miss_ps, tbd_miss_an, tbd_miss_pn];
        % -----------------------------------------------------------------------

        %% COMPARISONS WITH REZNIK ET AL. 2015
        % Here, we rearrange the results by taking into account only the trials
        % where the comparison tone was presented at the same intensity as the
        % standard one (i.e., relative difference/ relative intensity value = 0).
        % An additional difference in this analysis is that here we calculate the "1st Sound Louder response". This is done in order to visualize the data in the same way as Reznik et al. did.
        % The rest of the analysis is similar to the analysis described above.

        % Get rid of the trials where the standard intensity (column 2) was NOT the
        % same as the comparison one (column 3)
        same_PS = results_nomiss_ps; % for passive suprathreshold
        same_PS(same_PS(:, 2) ~= same_PS(:, 3), :) = [];

        same_AS = results_nomiss_as; % for active suprathreshold
        same_AS(same_AS(:, 2) ~= same_AS(:, 3), :) = [];

        same_AN = results_nomiss_an; % for active nearthreshold
        same_AN(same_AN(:, 2) ~= same_AN(:, 3), :) = [];

        same_PN = results_nomiss_pn; % for passive nearthreshold
        same_PN(same_PN(:, 2) ~= same_PN(:, 3), :) = []; %

        % Create new matrices to store only the same-intensity-trials and the
        % information needed to compare with renzik's studt

        % For AN
        reznik_comp_an = zeros(1, 3); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns
        reznik_comp_an(1, 1) = sum(same_AN(:, 4) == 0); % first column is sum of the "First tone louder" responses
        reznik_comp_an(1, 2) = length(same_AN); % third column is total number of times this intensity has been presented
        reznik_comp_an(1, 3) = (reznik_comp_an(1, 1) ./ reznik_comp_an(1, 2)) * 100; %Percent 1st Louder responses
        percent1ndLouder_AN = reznik_comp_an(1, 3) % Percent 1st Louder responses

        % For PN
        reznik_comp_pn = zeros(1, 3); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns
        reznik_comp_pn(1, 1) = sum(same_PN(:, 4) == 0); % first column is sum of the "First tone louder" responses
        reznik_comp_pn(1, 2) = length(same_PN); % third column is total number of times this intensity has been presented
        reznik_comp_pn(1, 3) = (reznik_comp_pn(1, 1) ./ reznik_comp_pn(1, 2)) * 100; %Percent 1st Louder responses
        percent1ndLouder_PN = reznik_comp_pn(1, 3) %Percent 1st Louder responses

        % For AS
        reznik_comp_as = zeros(1, 3); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns
        reznik_comp_as(1, 1) = sum(same_AS(:, 4) == 0); % first column is sum of the "First tone louder" responses
        reznik_comp_as(1, 2) = length(same_AS); % third column is total number of times this intensity has been presented
        reznik_comp_as(1, 3) = (reznik_comp_as(1, 1) ./ reznik_comp_as(1, 2)) * 100; %Percent 1st Louder responses
        percent1ndLouder_AS = reznik_comp_as(1, 3) %Percent 1st Louder responses

        % For PS
        reznik_comp_ps = zeros(1, 3); % create a new results matrix with as many rows as unique intensities delivered and 3 empty columns
        reznik_comp_ps(1, 1) = sum(same_PS(:, 4) == 0); % first column is sum of the "First tone louder" responses
        reznik_comp_ps(1, 2) = length(same_PS); % third column is total number of times this intensity has been presented
        reznik_comp_ps(1, 3) = (reznik_comp_ps(1, 1) ./ reznik_comp_ps(1, 2)) * 100; %Percent 1st Louder responses
        percent1ndLouder_PS = reznik_comp_ps(1, 3) %Percent 1st Louder responses

        fprintf(anal_logfile, '\r\n%s', '  ');
        fprintf(anal_logfile, '\r\n%s', ' Summary comparisons with Reznik et al. (2015) ');
        fprintf(anal_logfile, '\r\n%s', '  ');

        fprintf(anal_logfile, '\r\n%s', ['AS - P(1st Sound Louder) ' num2str(percent1ndLouder_AS)]);
        fprintf(anal_logfile, '\r\n%s', ['PS - P(1st Sound Louder) ' num2str(percent1ndLouder_PS)]);
        fprintf(anal_logfile, '\r\n%s', ['AN - P(1st Sound Louder) ' num2str(percent1ndLouder_AN)]);
        fprintf(anal_logfile, '\r\n%s', ['PN - P(1st Sound Louder) ' num2str(percent1ndLouder_PN)]);
        fprintf(anal_logfile, '\r\n%s', '  ');

        % %% SAVE
        cd(analysis_path)
        cd(analysis_scripts);
    end

    cd(analysis_path)
    eval(['save ' 'Summary_Excl_trials_Reviewer1 number_of_trials_to_be_excluded number_of_trials_to_be_excluded_x number_of_trials_to_be_excluded_miss']);

    cd(experiment_path)

    fclose(anal_logfile);
    AS = [];
    PS = [];
    AN = [];
    PN = [];

    RezAS = [];
    RezPS = [];
    RezAN = [];
    RezPN = [];
    %% WHOLE SAMPLE ANALYSIS
    % Creates files per condition including the data of the whole sample.
    % These matrices are later introduced in R for the fitting and simulations
    % procedures.

    % In the array "subjects" only the final sample needs to be included
    % Participants 2, 19, 25 excluded. Final Sample N = 28
    subjects = [1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 20 21 22 23 24 26 27 28 29 30 31];

    for subject = 1:length(subjects)
        load([analysis_path '/' num2str(subjects(subject), '%02d') '_Discrimination_AS_relative']);
        load([analysis_path '/' num2str(subjects(subject), '%02d') '_Discrimination_PS_relative']);
        load([analysis_path '/' num2str(subjects(subject), '%02d') '_Discrimination_AN_relative']);
        load([analysis_path '/' num2str(subjects(subject), '%02d') '_Discrimination_PN_relative']);

        % Load Reznik comp data
        load([analysis_path '/' num2str(subjects(subject), '%02d') '_percent1ndLouder_AS']);
        load([analysis_path '/' num2str(subjects(subject), '%02d') '_percent1ndLouder_PS']);
        load([analysis_path '/' num2str(subjects(subject), '%02d') '_percent1ndLouder_AN']);
        load([analysis_path '/' num2str(subjects(subject), '%02d') '_percent1ndLouder_PN']);

        results_forfit_as(:, 6) = subjects(subject);
        results_forfit_ps(:, 6) = subjects(subject);
        results_forfit_an(:, 6) = subjects(subject);
        results_forfit_pn(:, 6) = subjects(subject);

        AS = [AS; results_forfit_as];
        PS = [PS; results_forfit_ps];
        AN = [AN; results_forfit_an];
        PN = [PN; results_forfit_pn];

        percent1ndLouder_AS(:, 2) = subjects(subject);
        percent1ndLouder_PS(:, 2) = subjects(subject);
        percent1ndLouder_AN(:, 2) = subjects(subject);
        percent1ndLouder_PN(:, 2) = subjects(subject);

        RezAS = [RezAS; percent1ndLouder_AS];
        RezPS = [RezPS; percent1ndLouder_PS];
        RezAN = [RezAN; percent1ndLouder_AN];
        RezPN = [RezPN; percent1ndLouder_PN];
    end

    cd(analysis_path)
    eval(['save ' 'AS AS']);
    eval(['save ' 'PS PS']);
    eval(['save ' 'AN AN']);
    eval(['save ' 'PN PN']);
    eval(['save ' 'RezAS RezAS']);
    eval(['save ' 'RezPS RezPS']);
    eval(['save ' 'RezAN RezAN']);
    eval(['save ' 'RezPN RezPN']);
end
