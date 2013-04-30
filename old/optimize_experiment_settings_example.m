% Find best experimental setting (#lanes, #pools etc.)
% for a re-sequencing experiment example
% function tt = optimize_experiment_settings_example(abc)
% clear
auxStruct = struct;
auxStruct.maxNumTrials =  500;
auxStruct.sigma_pooling = 0.05; % variation in #reads between individuals in a pool
auxStruct.read_error_prob = 0.01; % prob. base read is different from true base
auxStruct.minFractionOfFailure = 0.05;
auxStruct.forceBBflag = 0;
auxStruct.single_region_length = 300; % average length of a region
auxStruct.read_length = 100; % length of each read
auxStruct.carrier_freq_vec = 1 ./ 100; % [100 500 1000]; % try different carrier allele frequencies
auxStruct.num_people_vec = [500 1000 2000 4000];
auxStruct.num_pools_vec = 10; % number of available pools/lanes
auxStruct.num_regions_vec = 1; % [1 10 100]; % determine regions size (and coverage)
auxStruct.num_people_in_pool_vec = 25; % [25 100 0.5]; % maximal number of people in one pool (fraction means frequency of people in pool)
auxStruct.num_barcodes = 1; % no barcodes
auxStruct.reads = 4*10^6 * auxStruct.read_length; % number used in the paper (here reads is measured in base-pairs)
auxStruct.simulate_direction = 'lanes_to_people';

auxStruct.known_str = str2word('_', auxStruct.simulate_direction, 1); % get the independent variable
auxStruct.unknown_str = str2word('_', auxStruct.simulate_direction, 3);
switch auxStruct.simulate_direction
    case 'lanes_to_people'
        auxStruct.known_var_vec = auxStruct.num_pools_vec;
    case 'people_to_lane'
        auxStruct.known_var_vec = auxStruct.num_people_vec;
end
N_opt = zeros(length(auxStruct.known_var_vec), length(auxStruct.num_regions_vec)); % optimal number of lanes/individuals
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
                cur_output_file = fullfile(['best_pooling_' num2str(auxStruct.num_regions_vec(j)) ...
                    '_regions_' num2str(auxStruct.known_var_vec(i)) '_' auxStruct.known_str '_' ...
                    num2str(auxStruct.freq) '_freq_' ...
                    num2str(auxStruct.num_people_in_pool) '_people_per_lane.mat']);
                auxStructFile = [cur_output_file(1:end-4) '_params.mat'];
                save(auxStructFile, 'auxStruct');
                N_opt(i,j) = num_lanes_to_max_N(auxStructFile, auxStruct.simulate_direction, ...
                    cur_output_file);
            end
        end
        output_lanes_file = fullfile(['estimated_N_opt_needed_freq_' ...
            num2str(auxStruct.freq) ...
            '_people_per_pool_' num2str(auxStruct.num_people_in_pool) '.txt']);
        save_experiment_settings(auxStruct, N_opt, output_lanes_file);
    end % loop on carrier frequency
end % loop on pool sparsity (num people per pool)