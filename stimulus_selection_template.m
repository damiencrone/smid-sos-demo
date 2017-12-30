% Example stimulus selection script for use with the Socio-Moral Image
% Database (SMID)
%   This script uses the SOS toolbox to draw sample(s) meeting a specific
%   set of criteria.
%
% Modfied:
%   December 2017 (Damien Crone)
% Requires:
%   The following toolboxes:
%       - SOS - see references below
%   The following files and folders:
%       - SMID_norms.csv - a .csv file containing SMID normative data
%       - SMID_img_properties.csv - a .csv file containing SMID image
%         proprty information
%       - a folder containing all stimuli
% For further information, see:
%   - The SOS toolbox overview paper:
%       Armstrong, B. C., Watson, C. E., & Plaut, D. C. (2012). SOS! An
%       algorithm and software for the stochastic optimization of stimuli.
%       Behavior Research Methods, 44(3), 675?705.
%       http://doi.org/10.3758/s13428-011-0182-9
%   - The SOS toolbox's online reference manual:
%       http://sos.cnbc.cmu.edu/commands.html
%   - The SMID OSF page:
%       https://osf.io/2rqad/
%
% This script was written and tested using MATLAB 2014b.


%% Preliminaries
clear all;
close all;
clc;


%% Specify paths

DIRS.stim_norms = 'SMID_norms.csv';
DIRS.stim_img_properties = 'SMID_img_properties.csv';
DIRS.stim = 'img/'; % Folder to retrieve images from
DIRS.out = ['output_', datestr(now, 'yyyy-mm-dd'), '/'];

% Make output directory and open log file
mkdir(char(DIRS.out));
diary([DIRS.out 'log.txt'])


%% Enter POPULATION specifications here

% Define image height to width ratio to limits and remove those outside
% (most images in the SMID lie between 0.6 and 0.8)
hwr_min = 0.6; % Lowest allowable HWR
hwr_max = 0.8; % Highest allowable HWR

% List any images to be excluded a priori
imgs_to_exclude = {...
    'b15_p385_17', ... image has multiple panels
    };
% Note that the this script should still execute properly if this cell
% array is empty


%% Enter SAMPLE  and OPTIMISATION specifications here

% Set size of each sample
set_size = 30;

% Maximium number of iterations for SOS
max_it = 2e6;

% Specify whether to exclude plots for specific variables
suppress_rating_freq_plots = true;


%% Load raw data and exclude images that are a priori ineligible

% Load raw data
SMID_NORMS = readtable(DIRS.stim_norms, 'ReadRowNames', true);
SMID_PROP = readtable(DIRS.stim_img_properties, 'ReadRowNames', true);
% Note that this table should only contain numeric variables (aside from
% the image IDs, which are assumed to be in the first column)

% Ensure tables are identically ordered
SMID_PROP = SMID_PROP(SMID_NORMS.Properties.RowNames, :);

% Merge tables
SMID = [SMID_NORMS.Properties.RowNames, SMID_NORMS, SMID_PROP];
SMID.Properties.VariableNames(1) = {'img'};

% Loop through and exclude specified images
n_before_apriori_exclusions = height(SMID);
for i = 1:length(imgs_to_exclude);
    str = imgs_to_exclude{i};
    excl = strcmpi(str, SMID.img);
    SMID(excl, :) = [];
end
disp(['Images excluded a priori = ', num2str(n_before_apriori_exclusions - height(SMID))]);

% Remove images with height to width ratios outside boundaries
excl = SMID.img_prop_height_width_ratio < hwr_min | ...
    SMID.img_prop_height_width_ratio > hwr_max;
SMID(excl, :) = [];
disp(['Images excluded based on HWR = ', num2str(sum(excl))]);


%% Save raw stimulus data as SOS-compatible text file

% Code first column as SOS string variable
col_names = SMID.Properties.VariableNames;
col_names{1} = [col_names{1}, '|s'];

% Code all remaining columns as SOS numeric variable
for i = 2:length(col_names);
    col_names{i} = [col_names{i}, '|f'];
end

% Convert data to cell array
C = table2cell(SMID);

% Initialise file
fileID = fopen('sos_data.txt', 'w') ;
formatSpec = [repmat('%s\t',1,length(col_names)-1), '%s\n'];

% Write re-formatted column names
fprintf(fileID, formatSpec, col_names{:});

% Modify formatting
formatSpec = ['%s\t', repmat('%5.3f\t',1,length(col_names)-1), '%5.3f\n'];

% Write rows of data
[nrows, ncols] = size(C);
for row = 1:nrows
    fprintf(fileID, formatSpec, C{row,:});
end

fclose(fileID);


%% Begin SOS

% Create a new optimizer
MoralSOS = sos('maxIt', max_it, 'reportInterval', 1000, 'stopFreezeIt', 1000, ...
    'statInterval', 5000, 'blockSize', 100);

% Load population
StimPop = population('sos_data.txt', ...
    'name','StimulusPopulation', ...
    'isHeader',true, ...
    'isFormatting',true ...
    );

% Define sample(s) - each element of the structure S corresponds to a
% different sample, which will produce a different text file output
S.moral   = sample(set_size, 'name', 'MoralSample');
S.immoral = sample(set_size, 'name', 'ImmoralSample');
S.neutral = sample(set_size, 'name', 'NeutralSample');

% Retrieve sample names
sample_names = fieldnames(S);

% Loop through samples to specify some additional settings
for i = 1:numel(sample_names)
    
    sn = sample_names{i};
    
    % Assign output file name
    S.(sn).outFile = [DIRS.out, S.(sn).name, '.txt'];
    
    % Assign stimulus population
    S.(sn).setPop(StimPop);
    
    % Add samples to optimizer
    MoralSOS.addSample(S.(sn));
    
end

% Fill initial samples with randomly selected stimuli
MoralSOS.initFillSamples();


%% Define constraints
%
% Here, each set of constraints is specified in separate tables of sample /
% variable / (etc.) combinations with one row per constraint. These tables
% are then used in the next code chunk to add the constraints to the
% optimizer.
%
% In this example, three sets of constraints are used:
%
%   (1) Hard constraints (defined in HARD_BOUND) that set upper and lower
%   limits on the moral ratings of images in the immoral and moral
%   conditions repectively
%
%   (2) Soft constraints (defined in SOFT_MATCH_1) that try and match
%   values in a given sample, on a given variable to a specific target. The
%   first row, for example, can be read as follows:
%
%       In the 'moral' sample, select images so that the 'moral_mean'
%       'mean' will equal 5
%
%   (3) Soft constraints (defined in SOFT_MATCH_2) that try and match two
%   sets of sample values to each other (i.e., to minimise the distance
%   between the two samples on those to variables).
%
%   (4) Soft constraints (defined in SOFT_DIST) that try to select images
%   so that in both conditions, the five moral foundation relevance
%   variables will be uniformly distributed, ensuring that the moral and
%   immoral samples contain diverse moral content.
%
% Other kinds of constraints can easily be added using a similar approach,
% for example matching correlations between variables within or across
% samples, minimizing the correlation between variables (e.g., between some
% moral variable and a potential confound like arousal), and so-on.

% (1) Define hard ceilings / floors on values for specific samples on
% specific variables
HARD_BOUND = readtable('constraint_tables/hard_bound.csv', 'ReadVariableNames', true);

% (2) Define soft constraints to match individual variables to target
% values
SOFT_MATCH_1 = readtable('constraint_tables/soft_match_1.csv', 'ReadVariableNames', true);

% (3) Define soft constraints to match variables to across samples
SOFT_MATCH_2 = readtable('constraint_tables/soft_match_2.csv', 'ReadVariableNames', true);

% (4) Define soft distributional constraints (i.e., try to return a
% solution where each of these variables are uniformly distributed)
SOFT_DIST = readtable('constraint_tables/soft_dist.csv', 'ReadVariableNames', true);


%% Set constraints

% This code chunk will loop through the tables defined above and add the
% specified constraints to the optimizer

% (1) Loop through hard ceiling / floor constraints and add each one to the
% optimizer
for i = 1:height(HARD_BOUND)
    
    % Construct constraint name
    constr_name = [ ...
        'hard_', char(HARD_BOUND.sample(i)), ...
        '_set_', char(HARD_BOUND.varname(i)), '_', char(SOFT_MATCH_1.stat(i)), ...
        '_to_', num2str(SOFT_MATCH_1.target(i)) ...
        ];
    
    MoralSOS.addConstraint(...
        'sosObj', MoralSOS, 'name', constr_name, ... Object and constraint names
        'constraintType', 'hard', 'fnc', char(HARD_BOUND.fnc(i)), ... Constraint type
        'sample1', S.(char(HARD_BOUND.sample(i))), ... Sample details
        's1ColName', char(HARD_BOUND.varname(i)), 'value', HARD_BOUND.value(i) ... Variable details
        );
    
end

% (2) Loop through soft match-to-value constraints and add each one to the
% optimizer
for i = 1:height(SOFT_MATCH_1)
    
    % Construct constraint name
    constr_name = [ ...
        'soft_', char(SOFT_MATCH_1.sample(i)), ...
        '_set_', char(SOFT_MATCH_1.varname(i)), '_', char(SOFT_MATCH_1.stat(i)), ... 
        '_to_', num2str(SOFT_MATCH_1.target(i)) ...
        ];
    
    MoralSOS.addConstraint(...
        ... Required options
        'sosObj', MoralSOS, 'name', constr_name, ... Object and constraint names
        'constraintType', 'soft', 'fnc', 'match1SampleVal', 'stat', SOFT_MATCH_1.stat(i), ... Constraint type
        'sample1', S.(char(SOFT_MATCH_1.sample(i))), ... Sample details
        's1ColName', char(SOFT_MATCH_1.varname(i)), 'targVal', SOFT_MATCH_1.target(i), ... Variable details
        ... Optional options
        'weight', SOFT_MATCH_1.weight(i), 'exponent', SOFT_MATCH_1.exp(i) ...
        );
    
end

% (3) Loop through soft match-value-across-samples constraints and add each
% one to the optimizer
for i = 1:height(SOFT_MATCH_2)
    
    % Construct constraint name
    constr_name = [ ...
        'soft_match_', char(SOFT_MATCH_2.sample1(i)), '_', char(SOFT_MATCH_2.varname1(i)), ...
        '_to_', char(SOFT_MATCH_2.sample2(i)), '_', char(SOFT_MATCH_2.varname2(i)), ...
        '_on_', char(SOFT_MATCH_2.stat(i)) ...
        ];
    
    MoralSOS.addConstraint(...
        ... Required options
        'sosObj', MoralSOS, 'name', constr_name, ... Object and constraint names
        'constraintType', 'soft', 'fnc', 'min', 'stat', SOFT_MATCH_2.stat(i), ... Constraint type
        'sample1', S.(char(SOFT_MATCH_2.sample1(i))), 'sample2', S.(char(SOFT_MATCH_2.sample2(i))), ... Sample details
        's1ColName', char(SOFT_MATCH_2.varname1(i)), 's2ColName', char(SOFT_MATCH_2.varname2(i)), ... Variable details
        'paired', false, ... Use pairwse or group-wise distances?
        ... Optional options
        'weight', SOFT_MATCH_2.weight(i), 'exponent', SOFT_MATCH_2.exp(i) ...
        );
    
end

% (4) Loop through soft distributional constraints and add each one to the
% optimizer
for i = 1:height(SOFT_DIST)
    
    % Construct constraint name
    constr_name = [ ...
        'soft_', char(SOFT_DIST.sample(i)), '_', char(SOFT_DIST.varname(i)), ...
        '_uniform_dist' ...
        ];
    
    
    MoralSOS.addConstraint(...
        ... Required options
        'sosObj', MoralSOS, 'name', constr_name, ... Object and constraint names
        'constraintType', 'soft', 'fnc', 'maxEnt', ... Constraint type
        'sample1', S.(char(SOFT_DIST.sample(i))), ... Sample details
        's1ColName', char(SOFT_DIST.varname(i)), 'pdSpread', 'allItems', ... Variable details
        ... Optional options
        'weight', SOFT_DIST.weight(i), 'nbin', 10, 'exponent', SOFT_DIST.exp(i) ...
        );
    
end


%% Run optimizer

% Set anneal schedule
MoralSOS.setAnnealSchedule('schedule', 'exp', 'blockSize',1000, ...
    'pDecrease', 0.46, 'pval',0.05);

% Initialise optimization
MoralSOS.optimize();


%% Handle output

% Loop through i samples
for i = 1:length(sample_names);
    
    % Create index variable for ith sample
    ind = char(sample_names(i));
    mkdir(char(strcat(DIRS.out, ind, '/')))
    
    % Write data
    S.(ind).writeData()
    
    % List images
    sample_imgs = S.(ind).data{1};
    
    % Loop through images
    for j = 1:length(sample_imgs);
        
        img_j = sample_imgs(j);
        img_fn = dir(char(strcat(DIRS.stim, img_j, '.*')));
        
        copyfile(char(strcat(DIRS.stim, img_fn.name)), ...
            char(strcat(DIRS.out, ind, '/')))
        
    end
    
end

% Save workspace
save([char(DIRS.out), 'workspace.mat']);

diary off;


%% Plot distributions for all variables

% Create plot output directories
mkdir(char(strcat(DIRS.out, 'plots/fig/')));
mkdir(char(strcat(DIRS.out, 'plots/jpg/')));

% Generate color matrix for samples
col_mat = hsv(numel(sample_names));

% Loop through unique variables
for var_num = 2:size(SMID, 2)
    
    % As a default, plot this variable
    is_var_to_plot = true;
    
    % Get variable name
    current_var = SMID.Properties.VariableNames(var_num);
    
    % Determine if variable is to be plotted
    current_var_char = char(current_var);
    nchar = numel(current_var_char);
    is_rating_freq = strcmp(current_var_char((nchar-1):end), '_n');
    
    if suppress_rating_freq_plots && is_rating_freq
        is_var_to_plot = false;
    end
    
    % If plot for variable is not supposed to be suppressed, go ahead and
    % plot it
    if is_var_to_plot
        
        % Add density plot for entire database
        [f,xi] = ksdensity(SMID.(char(current_var)));
        h = plot(xi, f, '--', 'color', [0.5, 0.5, 0.5]);
        
        hold all;
        
        % Loop through samples
        for sample_num = 1:numel(sample_names)
            
            current_sample = sample_names(sample_num);
            
            samp_col = col_mat(sample_num, :);
            
            % Get the first column of data (which contains the stimulus names)
            samp_stim = S.(char(current_sample)).data{1};
            
            % Get variable index
            var_names = [S.(char(current_sample)).header{:}];
            var_ind = find(strcmp(current_var, var_names));
            
            % Get values for current sample
            var_values = S.(char(current_sample)).data{var_ind};
            
            % Add density plot for current sample
            [f,xi] = ksdensity(var_values);
            h = [h, plot(xi, f, 'color', samp_col)];
            
            hold all;
            
        end % for sample_num
        
        % Add figure axis labels and legend
        xlabel(strrep(current_var, '_', ' '));
        ylabel('Density');
        legend(h, {'all', sample_names{:}});
        
        hold off
        
        % Save figure
        saveas(gcf, char(strcat(DIRS.out, 'plots/fig/', current_var)), 'fig');
        saveas(gcf, char(strcat(DIRS.out, 'plots/jpg/', current_var)), 'jpg');
        
    end % if is var to plot
    
end