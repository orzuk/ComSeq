% Compute the maximal sample size we can reconstruct for a given number of lanes
% This is done by performing a binary search on parameters of interest
%
% Input:
% auxStruct - structure with several experimental parameters
% simulate_direction - 'lanes_to_people' or 'people to lanes'
%
% Output:
% N_max - maximal number of individuals / minimal number of lanes allowing reliable reconstruction
% pooling_matrix - pooling matrix achieving reliable reconstruction
%
function [N_max pooling_matrix iters_performed] = ...
    num_lanes_to_max_N(auxStruct, simulate_direction, output_file, varargin)

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


