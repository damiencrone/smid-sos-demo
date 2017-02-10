% Example stimulus selection script for use with the Socio-Moral Image
% Database
%   This script uses the SOS toolbox to draw sample(s) meeting a specific
%   set of criteria.
%
% Modfied:
%   February 2017 (Damien Crone)
% Requires:
%   The following toolboxes:
%       - SOS - see references below
%   The following files and folders:
%       - stim_data.csv - a .csv file containing SMID normative data
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


%% Preliminaries
clear all; clc;


%% Specify paths

DIRS.stim_data = 'sos_data_unformatted.csv'; % Raw data
DIRS.stim = 'img/'; % Folder to retrieve images from
DIRS.out = ['output_', datestr(now, 'yyyy-mm-dd'), '/'];


%% Enter POPULATION specifications here

% Define image height to width ratio to limits and remove those outside
% (most images in the SMID lie between 0.6 and 0.8)
hwr_min = 0.6; % Lowest allowable HWR
hwr_max = 0.8; % Highest allowable HWR

% List any images to be excluded a priori
imgs_to_exclude = {...
    'b2_p17_2.jpg', ... contains famous people
    };
% Note that the this script should still execute properly if this cell
% array is empty


%% Enter SAMPLE  and OPTIMISATION specifications here

% Set size of each sample
set_size = 30;

% Maximium number of iterations for SOS
max_it = 2e6;


%% Load raw data and exclude images that are a priori ineligible

% Load raw data
T = readtable(DIRS.stim_data);
% Note that this table should only contain numeric variables (aside from
% the image IDs, which are assumed to be in the first column)

% Remove images with height to width ratios outside boundaries
excl = T.im_height_width_ratio < hwr_min | T.im_height_width_ratio > hwr_max;
T(excl, :) = [];

% Loop through and exclude specified images
for i = 1:length(imgs_to_exclude);
    str = imgs_to_exclude{i};
    excl = strcmpi(str, T.img);
    T(excl, :) = [];
end


%% Save raw stimulus data as SOS-compatible text file

col_names = T.Properties.VariableNames;
col_names{1} = [col_names{1}, '|s']; % Code as SOS string variable

for i = 2:length(col_names);
    col_names{i} = [col_names{i}, '|f']; % Code as SOS numeric variable
end

% Convert data to cell array
C = table2cell(T);

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
    'name','StimulusPopulation','isHeader',true,'isFormatting',true);

% Define sample(s) - each element of the structure S corresponds to a
% different sample, which will produce a different text file output
S.moral   = sample(set_size, 'name', 'MoralSample');
S.immoral = sample(set_size, 'name', 'ImmoralSample');

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

% Fill initial samples with random samples
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
%       In the 'moral' sample, select images so that the 'moral_all_avg'
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
%
% If using a large number of constraints, it may be easier to create these
% tables outside of matlab (e.g., in Excel), save them as a .csv file, and
% then load the table using the readtable() function.

% (1) Define hard ceilings / floors on values for specific samples on
% specific variables
HARD_BOUND = cell2table({...
    
      'moral', 'moral_all_avg',   'floor',       3.5;
    'immoral', 'moral_all_avg', 'ceiling',       2.5;

}, 'VariableNames', { ...

     'sample',       'varname',     'fnc', 'value'

});

% (2) Define soft constraints to match individual variables to target
% values
SOFT_MATCH_1 = cell2table({...
    
      'moral', 'moral_all_avg', 'mean',        5,        2,    5;
    'immoral', 'moral_all_avg', 'mean',        1,        2,    5;
      'moral',  'moral_all_sd', 'mean',        0,        1,    2;
    'immoral',  'moral_all_sd', 'mean',        0,        1,    2;

}, 'VariableNames', { ...

     'sample',       'varname', 'stat', 'target', 'weight', 'exp'

});

% (3) Define soft constraints to match variables to across samples
SOFT_MATCH_2 = cell2table({...
    
      'moral', 'immoral',      'harm_all_avg',      'harm_all_avg', 'mean',        2,    5;
      'moral', 'immoral',  'fairness_all_avg',  'fairness_all_avg', 'mean',        2,    5;
      'moral', 'immoral',   'ingroup_all_avg',   'ingroup_all_avg', 'mean',        2,    5;
      'moral', 'immoral', 'authority_all_avg', 'authority_all_avg', 'mean',        2,    5;
      'moral', 'immoral',    'purity_all_avg',    'purity_all_avg', 'mean',        2,    5;

}, 'VariableNames', { ...

    'sample1', 'sample2',          'varname1',          'varname2', 'stat', 'weight', 'exp'

});


% (4) Define soft distributional constraints (i.e., try to return a
% solution where each of these variables are uniformly distributed)
SOFT_DIST = cell2table({...
    
      'moral',      'harm_all_avg',        1,    2;
    'immoral',      'harm_all_avg',        1,    2;
      'moral',  'fairness_all_avg',        1,    2;
    'immoral',  'fairness_all_avg',        1,    2;
      'moral',   'ingroup_all_avg',        1,    2;
    'immoral',   'ingroup_all_avg',        1,    2;
      'moral', 'authority_all_avg',        1,    2;
    'immoral', 'authority_all_avg',        1,    2;
      'moral',    'purity_all_avg',        1,    2;
    'immoral',    'purity_all_avg',        1,    2;
    
}, 'VariableNames', { ...

     'sample',           'varname', 'weight', 'exp'
    
});


%% Set constraints

% This code chunk will loop through the tables defined above and add the
% specified constraints to the optimizer

% (1) Loop through hard ceiling / floor constraints and add each one to the
% optimizer
for i = 1:height(HARD_BOUND)
    
    % Construct constraint name
    constr_name = ['hard_', char(HARD_BOUND.sample(i)), '_set_', char(HARD_BOUND.varname(i)), '_', char(SOFT_MATCH_1.stat(i)), '_to_', num2str(SOFT_MATCH_1.target(i))];
    
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
    constr_name = ['soft_', char(SOFT_MATCH_1.sample(i)), '_set_', char(SOFT_MATCH_1.varname(i)), '_', char(SOFT_MATCH_1.stat(i)), '_to_', num2str(SOFT_MATCH_1.target(i))];
    
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
    constr_name = ['soft_match_', char(SOFT_MATCH_2.sample1(i)), '_', char(SOFT_MATCH_2.varname1(i)), '_', char(SOFT_MATCH_2.sample2(i)), '_', char(SOFT_MATCH_2.varname2(i)), '_', char(SOFT_MATCH_2.stat(i)),];
    
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
    constr_name = ['soft_', char(SOFT_DIST.sample(i)), '_', char(SOFT_DIST.varname(i)), '_uniform_dist'];
    

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

% Write data to text file
mkdir(char(DIRS.out));

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
        
        copyfile(char(strcat(DIRS.stim, '/', img_j)), ...
            char(strcat(DIRS.out, ind, '/')))
        
    end
    
end

% Save workspace
save([char(DIRS.out), 'workspace.mat'])