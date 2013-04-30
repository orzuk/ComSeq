% Prints the parameters used for simulations
function displaySimulationParameters(n_samples,k_pools,s,mean_reads,sigma_pipette,read_error_prob,varargin)

disp('simulation parameters: ')
disp(['number of individuals: ',num2str(n_samples)]);
disp(['# of pools: ',num2str(k_pools)]);
disp(['# carriers: ',num2str(s)]);
disp(['mean number of reads per SNP (see code for the interpretation): ',num2str(mean_reads)])
disp(['s.t.d of pipette error: ',num2str(sigma_pipette)]);
disp(['% read error: ',num2str(read_error_prob)])

for i=1:length(varargin)
  disp(['Remark: ',varargin{i}])
end
