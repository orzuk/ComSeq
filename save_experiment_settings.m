% Save in .txt the details of a pooling experiment
%
% Input:
% auxStruct - structure containing all experiment parameters
% N_opt - optimal number of lanes/people achieved
% K_opt - optimal ratio of lanes/pools 
% output_lanes_file - where to save results
%
% Output:
% R - structure with all experimental parameters
%
function R = save_experiment_settings(auxStruct, N_opt, K_opt, output_lanes_file)

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
    num2str(100 - 100*auxStruct.maxFractionOfFailure) '% of runs have no errors at all'];
R{5,1} = 'Parameters Used:';
R{6,1} = ['Region length = ' num2str(auxStruct.exon_length) ' basepairs'];
R{7,1} = ['Total Coverage Per Lane = ' num2str(auxStruct.reads) ' basepairs (aligned)'];
R{8,1} = ['Read Error = ' num2str(auxStruct.read_error_prob*100) '%'];
R{9,1} = ['Pooling Variation = ' num2str(auxStruct.sigma_pooling*100) '% (std. dev.)'];
R{10,1} = ['People per pool = ' num2str(auxStruct.num_people_in_pool)];
R{11,1} = ['Carrier Freq. = ' num2str(auxStruct.freq*100) '%'];
R{12,1} = ['Num. Barcodes = ' num2str(auxStruct.num_barcodes)];
R{13,1} = ['Sequencing cost per-lane = ' num2str(auxStruct.lane_cost) '$'];
R{14,1} = ['Amplification cost per-region = ' num2str(auxStruct.amplification_cost) '$'];
R{15,1} = '';
switch auxStruct.unknown_str
    case 'people'
        max_str = 'Maximal';
    case 'lanes'
        max_str = 'Minimal';
end
R{16,1} = [max_str ' Number of ' auxStruct.unknown_str ' Table:'];
R{17,1} = '';

offset = 18; % where to start putting the data
save_matrix = 0;
if(save_matrix) % old version: save matrix
    R{offset,1} = ['# Regions \ # ' auxStruct.known_str ' | '];
    for i=1:length(auxStruct.known_var_vec)
        R{offset,i+1} = num2str(auxStruct.known_var_vec(i));
    end
    R{offset+1,1} = '-----------------------------------------------------------------------------------------------------------------------------------';
    R{offset+1,2} = '---------------------------------------------------';
    for j=1:length(auxStruct.num_exons_vec)
        R{offset+1+j,1} = num2str(auxStruct.num_exons_vec(j));
        R{offset+1+j,1} = [R{offset+1+j,1} repmat(' ', 1, 20-length(R{offset+1+j,1})) '| '];
        %         R{15+j,1}(10) = '|'; % set at the same place
        %         R{15+j,1}(21) = ' '; % set at the same place
        for i=1:length(auxStruct.known_var_vec)
            R{offset+1+j,i+1} = num2str(N_opt(i,j));
        end
    end
else % new version: save vector of all (lanes, regions) options
    R{offset,1} = '# Regions |';
    R{offset,2} = ['# ' auxStruct.known_str ' |'];
    R{offset,3} = ['# ' auxStruct.unknown_str ' |'];
    R{offset,4} = '# Pools |'; %new! add also # pools 
    R{offset,5} = 'K |'; %new! add also K - the ratio between # pools and # lanes  
    R{offset,6} = 'Cost ($) |';
    R{offset,7} = 'Naive ($) |';
    R{offset,8} = 'Cost-Ratio |';
    R{offset,9} = 'Pools-Ratio';
    R{offset+1,1} = '-------------------------------------------------------------------------------------------------------------------------------------';
    ctr=1;
    for j=1:length(auxStruct.num_exons_vec)
        for i=1:length(auxStruct.known_var_vec)
            R{offset+1+ctr,1} = num2str(auxStruct.num_exons_vec(j));
            R{offset+1+ctr,1} = [R{offset+1+ctr,1}  repmat(' ', 1, 10-length(R{offset+1+ctr,1})) '| '];
            R{offset+1+ctr,2} = num2str(auxStruct.known_var_vec(i));
            R{offset+1+ctr,2} = [R{offset+1+ctr,2}  repmat(' ', 1, 9-length(R{offset+1+ctr,2})) '| '];
            R{offset+1+ctr,3} = num2str(N_opt(i,j));
            R{offset+1+ctr,3} = [R{offset+1+ctr,3}  repmat(' ', 1, 8-length(R{offset+1+ctr,3})) '| '];
            switch auxStruct.unknown_str
                case 'people'
                    num_lanes = auxStruct.known_var_vec(i); % num lanes is input 
                case 'lanes'
                    num_lanes = N_opt(i,j); % num lanes is output
            end
%             if(K_opt(i,j) == 0) 
%                 xxx = 999;
%             end
            num_pools = floor(num_lanes * auxStruct.num_barcodes / K_opt(i,j));
%             if( abs(num_pools -333) < 1)
%                 yyy = 2435;
%             end
            R{offset+1+ctr,4} = num2str(num_pools);
            R{offset+1+ctr,4} = [R{offset+1+ctr,4}  repmat(' ', 1, 8-length(R{offset+1+ctr,4})) '| '];
            
            R{offset+1+ctr,5} = num2str(K_opt(i,j)); % add K 
            R{offset+1+ctr,5} = [R{offset+1+ctr,5}  repmat(' ', 1, 2-length(R{offset+1+ctr,5})) '| '];

            auxStruct.region_length = auxStruct.num_exons_vec(j) * auxStruct.single_region_length;
            cur_cost = pooling_experiment_cost(auxStruct, num_lanes, num_pools); % New: compute cost (should be away from here)
            R{offset+1+ctr,6} = num2str(cur_cost);
            R{offset+1+ctr,6} = [R{offset+1+ctr,6}  repmat(' ', 1, 9-length(R{offset+1+ctr,6})) '| '];
%             if(cur_cost == 38885)
%                 xxx = 1234;
%             end
            naive_cost = pooling_experiment_cost(auxStruct, ...
                auxStruct.known_var_vec(i) / auxStruct.num_barcodes, auxStruct.known_var_vec(i)); % New: compute cost (should be away from here)% Now compute naive cost
            R{offset+1+ctr,7} = num2str(naive_cost);
            R{offset+1+ctr,7} = [R{offset+1+ctr,7}  repmat(' ', 1, 10-length(R{offset+1+ctr,7})) '| '];

            R{offset+1+ctr,8} = num2str(naive_cost / cur_cost);  % cost-ratio
            R{offset+1+ctr,8} = [R{offset+1+ctr,8}  repmat(' ', 1, 11-length(R{offset+1+ctr,8})) '| '];
            
            R{offset+1+ctr,9} = num2str(auxStruct.known_var_vec(i) / num_pools); % pools ratio 
            ctr=ctr+1;
        end
    end
end
bad_inds = find(cell2mat(str2nums_cell(R(offset+2:end,3))) == -1) + offset+1; % for these we don't yet have results 
good_inds = setdiff( 1:size(R,1), bad_inds);  R = R(good_inds,:); R{end+2,1} = ''; % create space 

my_mkdir(dir_from_file_name(output_lanes_file));
savecellfile(R, output_lanes_file);
