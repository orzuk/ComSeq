% This script computes the maximal sample size we can reconstruct for a given number of lanes, 
% or the minimal number of lanes sufficient to reconstruct a certain number of individuals.


% The struct auxStruct contains all relevant simulation parameters and then runs the program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auxStruct = struct;
auxStruct.simulate_direction = 'lanes_to_people'; % 'lanes_to_people' finds the maximal number of 
                                                  % individuals that can be reconstructed with zero
                                                  % errors for a given number of pools
                                                  
                                                  % 'people_to_lane' find the minimal number of lanes 
                                                  % needed to reconstruct a certain number of individuals

auxStruct.maxNumTrials =  500; % maximal number of simulations for each parameter setting
auxStruct.sigma_pooling = 0.05; % sigma of DNA preparation errors. 
auxStruct.read_error_prob = 0.01; % prob. base read is different from true base
auxStruct.minFractionOfFailure = 0.05; % maximal error level allowed. 
                                       %0.05 means that in at least 95% of the simulations one achieved 
                                       %zero error reconstruction
auxStruct.forceBBflag = 0; %if equals 'replace half by 2', then half of the
                           %non-zero entries in the genotype vector are replaces
                           %by the value 2, to simulate the "BB" case, .i.e homozygous rare-allele.
                           %Any other value for forceBBflag simulates the regular "AB" case.
                           %The default value of forceBBflag corresponds to the "AB" case.
auxStruct.single_region_length = 300; % average length of a region
auxStruct.read_length = 100; % length of each read
auxStruct.carrier_freq_vec = 1 ./ 100; %  % try different carrier allele frequencies
auxStruct.num_people_vec = [500 1000 2000 4000]; % try different numbers of individuals
auxStruct.num_pools_vec = 10; % number of available pools/lanes
auxStruct.num_regions_vec = 1; % [1 10 100]; % determine regions size (and coverage)
auxStruct.num_people_in_pool_vec = 25; % [25 100 0.5]; % maximal number of people in one pool (fraction means frequency of people in pool)
auxStruct.num_barcodes = 1; % "1" means no barcodes
auxStruct.reads = 4*10^6 * auxStruct.read_length; % number used in the paper (here reads are measured in base-pairs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% perform the simulation
optimize_experimental_setting(auxStruct)

% the output appears in a file - see directory.
