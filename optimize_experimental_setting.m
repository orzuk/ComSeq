% Computes the maximal sample size we can reconstruct for a given number of lanes, 
% or the minimal number of lanes sufficient to reconstruct a certain number of individuals.
% see optimize_experiment_setting_example.m for an example

function optimize_experimental_setting(auxStruct)
switch auxStruct.simulate_direction
    case 'lanes_to_people'
        auxStruct.known_var_vec = auxStruct.num_pools_vec;
    case 'people_to_lane'
        auxStruct.known_var_vec = auxStruct.num_people_vec;
end
N_opt = zeros(length(auxStruct.known_var_vec), length(auxStruct.num_regions_vec)); % optimal number of lanes/individuals

findFirstUnderscore = find(auxStruct.simulate_direction=='_');
auxStruct.known_str = auxStruct.simulate_direction(1:findFirstUnderscore(1)-1); % get the independent variable
auxStruct.unknown_str = auxStruct.simulate_direction(findFirstUnderscore(2)+1:end);


for p=1:length(auxStruct.num_people_in_pool_vec)  %loop on number of people per pool:
    auxStruct.num_people_in_pool = auxStruct.num_people_in_pool_vec(p);
    for c = 1:length(auxStruct.carrier_freq_vec) % loop on MAF
        auxStruct.freq = auxStruct.carrier_freq_vec(c);
        for i=1:length(auxStruct.known_var_vec) % loop on number of people. new: allow to loop over number of available lanes
            switch auxStruct.simulate_direction
                case 'lanes_to_people'
                    auxStruct.num_pools = auxStruct.known_var_vec(i);
                case 'people_to_lane'
                    auxStruct.num_people = auxStruct.known_var_vec(i);
            end
            if(auxStruct.num_people_in_pool < 1) % here take fraction of people in pool
                auxStruct.sqrtFlag = auxStruct.num_people_in_pool;
            else
                if(isfield(auxStruct, 'num_people'))
                    auxStruct.sqrtFlag = auxStruct.num_people_in_pool / auxStruct.num_people;
                else
                    auxStruct.sqrtFlag = auxStruct.num_people_in_pool;
                end
            end
            
            for j=1:length(auxStruct.num_regions_vec) % loop on region length covered
                auxStruct.region_length = auxStruct.num_regions_vec(j) * auxStruct.single_region_length;
                cur_output_file = ['best_pooling_' num2str(auxStruct.num_regions_vec(j)) ...
                    '_regions_' num2str(auxStruct.known_var_vec(i)) '_' auxStruct.known_str '_' ...
                    num2str(auxStruct.freq) '_freq_' ...
                    num2str(auxStruct.num_people_in_pool) '_people_per_lane.mat'];
                auxStructFile = [cur_output_file(1:end-4) '_params.mat'];
                save(auxStructFile, 'auxStruct');
                N_opt(i,j) = num_lanes_to_max_N(auxStructFile, auxStruct.simulate_direction, ...
                    cur_output_file);
            end
        end
        output_lanes_file = ['estimated_N_opt_needed_freq_' num2str(auxStruct.freq) ...
            '_people_per_pool_' num2str(auxStruct.num_people_in_pool) '.txt'];
        save_experiment_settings(auxStruct, N_opt, output_lanes_file);
    end % loop on carrier frequency
end % loop on pool sparsity (num people per pool)


function [N_max pooling_matrix iters_performed] = ...
    num_lanes_to_max_N(auxStruct, simulate_direction, output_file, varargin)

% this function performs the binary search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAX_LANES = 500; % maximal number currently reasonable

if(ischar(auxStruct)) % input is file name
    load(auxStruct);
end
if(isfield(auxStruct, 'maxNumTrials'))
    iters = auxStruct.maxNumTrials;
else
    iters = 1000;  % how many simulations for each parameters settings
end
if(isfield(auxStruct, 'minFractionOfFailure'))
    alpha = auxStruct.minFractionOfFailure;
else
    alpha = 0.05; % fraction of reconstructions we allow to fail
end

switch simulate_direction
    case 'lanes_to_people'
        N_min = auxStruct.num_pools; % naive approach: one lane per person
        N_max = min(10000, floor(auxStruct.num_pools / auxStruct.freq)); % best possible approach (depends on sparsity)
    case 'people_to_lanes'
        N_max = ceil(auxStruct.num_people / auxStruct.num_barcodes); % naive approach: one lane per person
        N_max = min(N_max, MAX_LANES); % don't allow too many lanes
        N_min = floor(auxStruct.freq * auxStruct.num_people * ...
            log(auxStruct.num_people) / auxStruct.num_barcodes); % best possible approach (logarithmic, depends on sparsity)
end
N_interval_searched = [N_min N_max];
mean_reads = auxStruct.reads / ...
    (auxStruct.region_length .* auxStruct.num_barcodes); % split reads into regions

while (N_min < N_max-1) % we can stop when difference is one
    N_mid = floor((N_min + N_max)/2)
    switch simulate_direction
        case 'lanes_to_people'
            num_pools = auxStruct.num_pools * auxStruct.num_barcodes; % consider barcoding
            num_people = N_mid;
        case 'people_to_lanes'
            num_pools = N_mid  * auxStruct.num_barcodes;
            num_people = auxStruct.num_people;
    end
    if(auxStruct.sqrtFlag > 1) % here we set a constant number of people in a pool
        auxStruct.UseSqrtFlag = auxStruct.sqrtFlag / num_people;
    else
        auxStruct.UseSqrtFlag = auxStruct.sqrtFlag;
    end
    
    num_reconstruction_errors = zeros(iters,1); % make sure this is set to zero again each time
    for i=1:iters
        if(mod(i,10) == 0)
            cur_iter = i
        end
        [x, fractionalOutput, discreteOutput ...
            num_reconstruction_errors(i) pooling_matrix] = ...
            simulateCSseq(num_people, num_pools, auxStruct.freq, mean_reads, ...
            auxStruct.sigma_pooling, auxStruct.read_error_prob, ...
            auxStruct.UseSqrtFlag, auxStruct.forceBBflag);
        current_failures = sum(num_reconstruction_errors > 0); % we may decide to abort online (if enough errors have passed)
        current_successes = i - current_failures;
        break_flag = decide_success_or_failure_online( ...
            current_successes, current_failures, alpha, iters);
        if(break_flag > 0)
            success_flag = break_flag -1;
            stopped_after_iters = i
            break;
        end
    end
    if(i == iters)
        success_flag = sum(num_reconstruction_errors > 0) < alpha * iters;
    end
    iters_performed = i;
    switch simulate_direction
        case 'people_to_lanes' % reverse direction (increase lanes)
            success_flag = 1-success_flag;
    end
    if(success_flag) % increase people/decrease lanes
        N_min = N_mid;
    else % decrease people/increase lanes
        N_max = N_mid;
    end
end
if(exist('output_file', 'var')) % save output to file
    save(output_file, 'auxStruct', 'N_max', 'pooling_matrix', 'N_interval_searched');
end

function break_flag=decide_success_or_failure_online( current_successes, current_failures, alpha, iters)
% Internal function for deciding when to stop trying (if number of failures
% is too large). Currently be conservative: stop if at least alpha*N
% failures (we can do better by getting confidence intervals on the
% observed alpha)

break_method = 'binomial_confidence_interval'; % if we're confidence number is out of range 
switch break_method
    case 'reach_alpha_failures'
        if(current_failures >= alpha*iters)
            break_flag = 1;
        else
            break_flag = 0;
        end
        
    case 'binomial_confidence_interval' % can be improved by lowering variance on the extremes
        N = current_successes + current_failures;
        p = current_failures / N; % estimator of frequency of errors
        beta = 0.001; % probability we're outside the confidence interval 
        confidence_sigma = norminv(1-beta/2) * sqrt(1/(4*N));
        confidence_interval = [p - confidence_sigma p + confidence_sigma]; 
        
        break_flag = 0; 
        if(alpha < confidence_interval(1)) % here we've failed 
            break_flag = 1;
        end
        if(alpha > confidence_interval(2)) % here we've succeeded
           break_flag = 2;  
        end
end

function R = save_experiment_settings(auxStruct, N_opt, output_lanes_file)
% Save in .txt the details of a pooling experiment

if(~isfield(auxStruct, 'num_exons_vec'))
    auxStruct.num_exons_vec = auxStruct.num_regions_vec;
end
if(~isfield(auxStruct, 'exon_length'))
    auxStruct.exon_length = auxStruct.single_region_length;
end
R = cell(16+length(auxStruct.num_exons_vec),1+length(auxStruct.num_people_vec));
R{1,1} = 'Optimal Parameters Settings. Number of Lanes Needed for Accurate Reconstruction';
R{2,1} = '--------------------------------------------------------------------------------';
R{3,1} = '';
R{4,1} = ['Reconstruction Success Criteria: ' ...
    num2str(100 - 100*auxStruct.minFractionOfFailure) '% of runs have no errors at all'];
R{5,1} = 'Parameters Used:';
R{6,1} = ['Region length = ' num2str(auxStruct.exon_length) ' basepairs'];
R{7,1} = ['Total Coverage Per Lane = ' num2str(auxStruct.reads) ' basepairs (aligned)'];
R{8,1} = ['Read Error = ' num2str(auxStruct.read_error_prob*100) '%'];
R{9,1} = ['Pooling Variation = ' num2str(auxStruct.sigma_pooling*100) '% (std. dev.)'];
R{10,1} = ['People per pool = ' num2str(auxStruct.num_people_in_pool)];
R{11,1} = ['Carrier Freq. = ' num2str(auxStruct.freq*100) '%'];
R{12,1} = '';
switch auxStruct.unknown_str
    case 'people'
        max_str = 'Maximal';
    case 'lanes'
        max_str = 'Minimal';
end
R{13,1} = [max_str ' Number of ' auxStruct.unknown_str ' Table:'];
R{14,1} = '';
R{15,1} = ['# Regions \ # ' auxStruct.known_str ' | '];
for i=1:length(auxStruct.known_var_vec)
    R{15,i+1} = num2str(auxStruct.known_var_vec(i));
end
R{16,1} = '--------------------------------------------------------------------------------';
for j=1:length(auxStruct.num_exons_vec)
    R{16+j,1} = num2str(auxStruct.num_exons_vec(j));
    R{16+j,1} = [R{16+j,1} repmat(' ', 1, 19-length(R{16+j,1})) '| '];
    %         R{15+j,1}(10) = '|'; % set at the same place
    %         R{15+j,1}(21) = ' '; % set at the same place
    for i=1:length(auxStruct.known_var_vec)
        R{16+j,i+1} = num2str(N_opt(i,j));
    end
end
%     save_expression_mat_as_text_generic(num2cell(num_people_vec), ...
%         num2cell(num_exons_vec)', N_opt', ...
%         output_lanes_file, 1); % need to transpose lanes and people !!!!
%
savecellfile(R, output_lanes_file);

function savecellfile(table, outFile, delimiter, save_mode, varargin)
% Save a cell array into a tab-delimited .txt file (Written by Assif Yitzhaky)
%
% Input:
% table - the cell array to be saved
% outfile - where to save the cell
% delimiter - what delimiter to use when saving to file (default is tab '\t')
% save_mode - whether to save strings without the '%s' (default, 0) or with
% (1). In the default state, '\t' will be converted into tabs, and '\n' to new lines
%

table = empty_cell_to_empty_str(table); % convert empty cells to empty strings 

if((nargin < 3) || isempty(delimiter))
    delimiter = '\t';
end
if(nargin < 4) 
    save_mode = 0;
end
if(~isnumeric(save_mode))
    save_mode= 0;
end
fout=fopen(outFile, 'wt');

% saving_lines = size(table,1)
for i=1:size(table, 1)
   % i_is = i
    if(mod(i, 1000) == 0)
        sprintf('save line %ld out of %ld',  i, size(table,1))
    end
    for j=1:size(table, 2)-1
        if isempty(table{i,j})
            fprintf(fout, delimiter);
            continue;
        elseif isnumeric(table{i,j})
            specifier=['%d' delimiter];
        else
            specifier=['%s' delimiter];
        end
        fprintf(fout, specifier, table{i,j});
    end
    if isnumeric(table{i, size(table, 2)})
        specifier='%d';
    else
        specifier='%s';
    end
    if(~isempty(strmatch('%s', specifier))  && (~save_mode) && isempty(strfind(table{i, size(table, 2)}, 'href=' ))) % we don't use the %s - to avoid tabs being lost .. - but then the problem is that \n are turned into new lines ..
        fprintf(fout, table{i, size(table, 2)} );
    else
        fprintf(fout, specifier, table{i, size(table, 2)} );
    end
    fprintf(fout, '\n');
end
fclose(fout);

function c_str = empty_cell_to_empty_str(c)
% Convert empty cells in a cell array to empty strings
c_str = c;
for i=1:size(c_str,1)
    for j=1:size(c_str,2)
        if(isempty(c_str{i,j}))
            c_str{i,j} = '';
        end
    end
end
